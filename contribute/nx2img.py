import networkx as nx
import sys
import os

STRUCT = 'structure'
TRIMMED_MCS_FRAG = 'trimmed_mcs_frag'
ORIGINAL_MCS_FRAG = 'original_mcs_frag'

def guess_format(fname):
    format = None
    root, ext = os.path.splitext(fname)
    if ext:
        format = ext[1:].lower() 
    
    return format


class TableRender:
    
    def __init__(self, table_border=0, table_cellborder=0, table_cellspacing=0):
        # FIXME
        self._rows = []
        self._table_border = table_border
        self._table_cellborder = table_cellborder
        self._table_cellspacing = table_cellspacing
    
    def _processTextCell(self, text):
        return '<TD>%s</TD>'%text
    
    def _processImageCell(self, text):
        return '<TD><IMG SRC="%s"></IMG></TD>'%text    
    
    def _processRow(self, cells):
        return '<TR>%s</TR>'%cells
    
    
    def _processTable(self):
        rows = ''.join(self._rows)
        return '<<TABLE BORDER="%d" CELLBORDER="%d" CELLSPACING="%d">%s</TABLE>>'%(self._table_border, 
                                                                                   self._table_cellborder,
                                                                                   self._table_cellspacing,
                                                                                   rows)
    
    def addText(self, text):
        if not isinstance(text, list):
            text = [text]
        cells = ''.join([ self._processTextCell(e) for e in text])
        row = self._processRow(cells)
        self._rows.append(row)
        
    
    def addImage(self, image):
        if not isinstance(image, list):
            image = [image]        
            
        cells = ''.join([self._processImageCell(e) for e in image])
        row = self._processRow(cells)
        self._rows.append(row)
    
    def getHtml(self):
        return self._processTable()
        
    
    
class DotRender:
    def __init__(self, G, img_generator, 
                 node_attributes, 
                 edge_attributes, 
                 save_image=True,
                 img_format='svg',
                 font_size=40,
                 ):
        
        self._G = G
        D = self.nx2dot(G)
        self._D = D
        self._img_list = []
        self._img_generator = img_generator
        self._node_attributes = node_attributes
        self._edge_attributes = edge_attributes
        self._save_image = save_image
        self._img_format = img_format
        self._basename = None
        self._font_size = font_size
    
    def nx2dot(self, G):
        H = G.copy()
        
        for node1_id, node2_id in H.edges():
            removed_list = []
            for key, value in H.get_edge_data(node1_id, node2_id).items():
                if isinstance(value, dict):
                    removed_list.append(key)
            for key in removed_list:
                del H.edge[node1_id][node2_id][key]
        
        D = nx.to_pydot(H)
        return D
    def run(self, filename, format=None):
        
        self._basename, suffix = os.path.splitext(os.path.basename(filename))
        
        self._img_list = []
        
        self._renderNodes()
        self._renderEdges()
   
        self._D.set_shape_files(self._img_list)
        self._D.write(filename, format=format)
        
        if not self._save_image:
            for image in self._img_list:
                os.remove(image)   
                
    def _getImgFname(self, name):
        img_fname = '%s_%s.%s'%(self._basename, 
                                name, 
                                self._img_format)
        return img_fname
    
    def _renderNodeStructure(self, dot_node, smiles, img_fname):
        self._img_generator.generate(smiles, img_fname)
    
    def _renderEdgeStructure(self, source, dest, attr_name, 
                             source_img_fname, dest_img_fname):
        source_id = self._getNodeName(source)
        dest_id = self._getNodeName(dest)
        source_mol = self._G.edge[source_id][dest_id][attr_name][source_id]           
        img_generator.mol2svg(source_mol, source_img_fname)
        dest_mol = self._G.edge[source_id][dest_id][attr_name][dest_id]
        img_generator.mol2svg(dest_mol, dest_img_fname)
            
        
    
    def _getNodeName(self, node):
        # networkx node id is a hex string
        # pydot seems to add '"' to id that begins with digit
        name = node.get_name()
        name = name.strip('"')
        return name
        
    def _renderNodes(self):
        for node in self._D.get_node_list():    
            
            node.set_fontsize(self._font_size)   
            name = self._getNodeName(node)
            node_index = node.get_sequence()
            render = TableRender()
            for attr_name  in self._node_attributes:
                attr = node.get(attr_name)
                if attr:
                    if attr_name == 'SMILES':
                        img_fname = self._getImgFname("node%d"%node_index)
                        
                        self._renderNodeStructure(node, attr, img_fname)
                        render.addImage(img_fname)
                        self._img_list.append(img_fname)  
                    elif attr_name == 'title':
                        render.addText('%s'%attr)
                    else:
                        render.addText('%s:%s'%(attr_name, attr))
        
            node.set_label(render.getHtml())        
    
    def _renderEdges(self):
        simi  = [float(edge.get('similarity')) for edge in self._D.get_edges()]
        scale = 1.0 / max( simi )           
        for i, edge in enumerate(self._D.get_edge_list()): 
            
            edge_index = i + 1
            edge.set_fontsize(self._font_size)        
            
            try :
                partial_ring = edge.get("partial_ring")
                if partial_ring is None:
                    partial_ring = 0
            except ValueError :
                partial_ring = 0
            saturation = float(edge.get("similarity")) * scale
            saturation = 0.0 if (saturation < 0) else (1.0 if (saturation > 1) else saturation)
            edge.set_color("0.8396,%f,0.8" % saturation)
            edge.set_weight(saturation)
            edge.set_penwidth(6.0)
            if saturation < 0.01 or partial_ring:
                edge.set_style("dashed")
            
            source_id = edge.get_source()
            dest_id = edge.get_destination()
            source = self._D.get_node(source_id)[0]
            dest = self._D.get_node(dest_id)[0]
            source_title = source.get('title')
            dest_title = dest.get('title')  
            
            render = TableRender(table_border=0)
            render.addText(["%s : %s"%(source_title, dest_title)])
            
            edge_index = edge.get_sequence()
            for attr_name in self._edge_attributes:
                try:
                    attr = edge.get(attr_name)
                    if attr is not None:
                        if attr_name == 'SMILES':

                            img_fname = self._getImgFname("edge%d"%(edge_index))
                            if edge.get("layout_mcs"):
                                attr = edge.get("layout_mcs")
                            self._img_generator.generate(attr, img_fname)
                            render.addImage(img_fname)
                            self._img_list.append(img_fname)    
                        else:
                            render.addText("%s:%s"%(attr_name, attr))
                        
                except KeyError, e:
                    print e
                    pass        
            
            edge.set_label(render.getHtml()) 
            
    


if __name__ == '__main__':
    usage = \
    """usage: %prog [options] <input.pkl> <output.svg>
    A script to generate an image file from a networkx graph pickle file. 
    Each node must have a SMILES attribute that represents the associated molecule.
    """
    from optparse import OptionParser
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--save",
                      action="store_true", dest="save", default=False,
                      help="Save intermediate images for individual molecule.")

    parser.add_option("-f", "--format",
                      action="store", dest="format", default='svg',
                      choices=["bmp", "jpg", "jpeg", "png", "ppm", "tiff", "xbm", "xpm", "svg"],
                      help="Intermediate image format.")
    
    parser.add_option("-n", "--node_attribute",
                      action="append", dest="node_attributes", default=[],
                      help="node attribues to be rendered. Default attributes incluing 'SMILES' and 'title'.")    
    parser.add_option("-e", "--edge_attribute",
                      action="append", dest="edge_attributes", default=[],
                      help="node attribues to be rendered. 'SMILES' is the default attribute.")        
    
    
    options, arg = parser.parse_args()
    if len(arg) != 2:
        parser.print_help()
        sys.exit(1)

    if len(options.node_attributes) == 0:
        options.node_attributes.append('SMILES')
        options.node_attributes.append('title')
        
    if len(options.edge_attributes) == 0:
        options.edge_attributes.append('SMILES')
  
        
    input_pickle = arg[0]
    output_fname = arg[1]
    
    format = guess_format(output_fname)
    if options.format == 'svg' and format == 'svg' and not options.save:
        print "WARNING:the output svg file depends on intermediate svg files, please use -s option."
        sys.exit(0)

    G = nx.read_gpickle(input_pickle)    
    
    
  
    try:
        from imggenerator import VirtualX
        from imggenerator import OpeneyeImgGenerator
        img_generator = OpeneyeImgGenerator()

        render = DotRender(G, img_generator, 
                       options.node_attributes,
                       options.edge_attributes,
                       save_image=options.save,
                       img_format=options.format,
                       )
    except ImportError, e :
        pass
    try:
        from imggenerator import VirtualX
        from imggenerator import SchrodImgGenerator
        img_generator = SchrodImgGenerator()
        render = DotRender(G, img_generator, 
                       options.node_attributes,
                       options.edge_attributes,
                       save_image=options.save,
                       img_format=options.format,
                       )
    except ImportError, e :
        pass

    print "rendering..."
    render.run(output_fname, format=format)
    

