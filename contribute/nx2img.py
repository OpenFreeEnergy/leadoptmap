import networkx as nx
import sys
import os
from imggenerator import SchrodImgGenerator, VirtualX

STRUCT = 'structure'
def guess_format(fname):
    format = None
    root, ext = os.path.splitext(fname)
    if ext:
        format = ext[1:].lower() 
    
    return format

def render_smiles(smiles, img_fname):
    img_generator.generate(smiles, img_fname)
    return '<TR><TD><IMG SRC="%s"/></TD></TR>'%img_fname

def render_structure(mol, img_fname):
    img_generator.mol2svg(mol, img_fname)
    return '<TR><TD><IMG SRC="%s"/></TD></TR>'%img_fname

def render_attr(attr_name, text):
    return render_text("%s: %s"%(attr_name, text))

def render_text(text):
    return "<TR><TD>%s</TD></TR>"%(text)

class DotRender:
    def __init__(self, G, img_generator, 
                 node_attributes, 
                 edge_attributes, 
                 save_image=True,
                 img_format='svg',
                 font_size=40,
                 align=False):
        
        self._G = G
        self._align = align
        if self._align:
            self._alignGraph()
        
        D = nx.to_pydot(G)
        self._D = D
        self._img_list = []
        self._img_generator = img_generator
        self._node_attributes = node_attributes
        self._edge_attributes = edge_attributes
        self._save_image = save_image
        self._img_format = img_format
        self._basename = None
        self._font_size = font_size
    
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
        if self._align:
            node_id = self._getNodeName(dot_node)
            mol = self._G.node[node_id][STRUCT]
            return render_structure(mol, img_fname)
        else:
            return render_smiles(smiles, img_fname)
    
    def _renderEdgeStructure(self, source, dest, smiles, img_fname):
        if self._align:
            source_id = self._getNodeName(source)
            dest_id = self._getNodeName(dest)
            mol = self._G.edge[source_id][dest_id][STRUCT]
            return render_structure(mol, img_fname)
        else:
            return render_smiles(smiles, img_fname)
        
    
    def _getNodeName(self, node):
        # networkx node id is a hex string
        # pydot seems to add '"' to id that begins with digit
        name = node.get_name()
        name = name.strip('"')
        return name
        
    def _renderNodes(self):
        for i, node in enumerate(self._D.get_node_list()):    
            node_index = i + 1
            node.set_fontsize(self._font_size)   
            name = self._getNodeName(node)
            render_string = ''
                       
            for attr_name  in self._node_attributes:
                attr = node.get(attr_name)
                if attr:
                    if attr_name == 'SMILES':
                        img_fname = self._getImgFname("node_%06d"%node_index) 
                        render_string = render_string + \
                            self._renderNodeStructure(node, attr, img_fname)
                        self._img_list.append(img_fname)                  
                    else:
                        render_string = render_string + render_attr(attr_name, attr)
                        
                
        
            node.set_label('<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0">%s</TABLE>>'%(render_string))        
    
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
            
            render_string = render_text("[%s,%s]"%(source_title, dest_title))            
            
            for attr_name in self._edge_attributes:
                try:
                    attr = edge.get(attr_name)
                    if attr:
                        if attr_name == 'SMILES':
                            img_fname = self._getImgFname("edge_%06d"%edge_index)
                            render_string = render_string + \
                                self._renderEdgeStructure(source, dest, attr, img_fname)
                            self._img_list.append(img_fname)                       
                        else:
                            render_string = render_string + render_attr(attr_name, attr)
                        
                except KeyError, e:
                    print e
                    pass        
            
            edge.set_label('<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0">%s</TABLE>>'%(render_string)) 
    
    def _alignGraph(self):
        
        from schrodinger.structutils.analyze import evaluate_smarts_canvas
        from schrodinger.infra import canvas2d
        from operator import itemgetter
        from schrodinger.application.canvas.base import ChmLicenseShared_isValid
        from schrodinger.application.canvas.utils import get_license
        if not ChmLicenseShared_isValid():
            canvas_license = get_license("LICENSE_SHARED")           
        
        IGNORE_HYDROGEN = 0
        RESCALE = 2
        FIXUP = True        
        
        def smiles2mol(smiles):
            mol = canvas2d.ChmMol.fromSMILES(smiles)    
            # need to generate 2D coordinates, otherwise 
            # the following function call to generateFromTemplateAndApply
            # will not have any effects
            canvas2d.Chm2DCoordGen.generateAndApply(mol, IGNORE_HYDROGEN)
            return mol
        
        def evaluate_smarts(mol2d, smarts):
            mol3d = canvas2d.convertChmMoltoSWIG(mol2d)
            return evaluate_smarts_canvas(mol3d, smarts, start_index=0)            


        aligned = set()
        
        # loop over connected component of the graph
        for node_list in nx.connected_components(self._G):

            # sort the nodes by their degrees.
            sorted_nodes = sorted(self._G.degree(node_list).items(), \
                                  key=itemgetter(1), reverse=True)
            node_list = [ n[0] for n in sorted_nodes]

            # select the node with biggest degree as reference node
            ref_node_id = sorted_nodes[0][0]
            ref_node = self._G.node[ref_node_id]
            ref_node[STRUCT] = smiles2mol(ref_node['SMILES'])
            aligned.add(ref_node_id)

            # doing breadth-first search to get the edge list
            edge_list = []
            for id1, id2 in nx.bfs_edges(self._G, ref_node_id):
                if id1 < id2:
                    edge_list.append((id1, id2))
                else:
                    edge_list.append((id2, id1))
            
            edge_set = set(edge_list)

            # adding unvisited edges to edge_list 
            for id1, id2 in self._G.edges(node_list):
                if id1 > id2:
                    tmp = id1
                    id1 = id2
                    id2 = tmp
                    
                if (id1, id2) in edge_set:
                    continue
                edge_list.append((id1, id2))
                
                
            # loop over edge list to align 2D structure. 
            for node1_id, node2_id in edge_list:
                node1 = self._G.node[node1_id]
                node2 = self._G.node[node2_id]
    
                if STRUCT not in node1:
                    node1[STRUCT] = smiles2mol(node1['SMILES'])
                if STRUCT not in node2:
                    node2[STRUCT] = smiles2mol(node2['SMILES'])         
                    
                mol1 = node1[STRUCT]
                mol2 = node2[STRUCT] 
                edge_data = self._G.get_edge_data(node1_id, node2_id)
                # when the backend is adjusted, we should have a pair
                # of SMARTS that will match one node of the edge
                # for now we just use the SMILES and assume it will match
                # both nodes.
                smarts = edge_data['SMILES'] 
                
                ma1 = evaluate_smarts(mol1, smarts)
                ma2 = evaluate_smarts(mol2, smarts)
                
                if len(ma1) > 0 and len(ma2) > 0:
                    if node1_id in aligned and node2_id in aligned:
                        pass
                    elif node1_id not in aligned:
                        # use node2 as reference
                        canvas2d.Chm2DCoordGen.generateFromTemplateAndApply(\
                            mol1, mol2,
                            ma1[0], ma2[0], 
                            IGNORE_HYDROGEN,
                            RESCALE,
                            FIXUP)
                        aligned.add(node2_id)
                    elif node2_id not in aligned:
                        # use node1 as reference
                        canvas2d.Chm2DCoordGen.generateFromTemplateAndApply(\
                            mol2, mol1,
                            ma2[0], ma1[0], 
                            IGNORE_HYDROGEN,
                            RESCALE,
                            FIXUP)  
                        aligned.add(node1_id)
                    else:
                        print "can not align (%s, %s)"%(node1_id, node2_id)
                else:
                    print "(%s, %s) mismatched"%(node1['title'], node2['title'])
                
                #align MCS
                template_mol = None
                template_core = None
                if len(ma1) > 0:
                    template_mol = mol1
                    template_core = ma1[0]
                elif len(ma2) > 0:
                    template_mol = mol2
                    template_core = ma2[0]   
                else:
                    print "unable to align MCS for edge(%s, %s)"%(node1_id, node2_id)
                
                mcs_frag = smiles2mol(edge_data['SMILES'])
                            
                if template_mol and len(template_core) > 0:
    
                    canvas2d.Chm2DCoordGen.generateFromTemplateAndApply(\
                        mcs_frag, template_mol,
                        range(len(template_core)), template_core, 
                        IGNORE_HYDROGEN,
                        RESCALE,
                        FIXUP)   
                    
                self._G.edge[node1_id][node2_id][STRUCT] = mcs_frag 
            
   



if __name__ == '__main__':
    usage = \
    """usage: %prog [options] <input.pkl> <output.png>
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
    parser.add_option("-a", "--align",
                      action="store_true", dest="align", default=False,
                      help="Align 2D image.")
    
    
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
    
    
  
    img_generator = SchrodImgGenerator()

    render = DotRender(G, img_generator, 
                       options.node_attributes,
                       options.edge_attributes,
                       save_image=options.save,
                       img_format=options.format,
                       align=options.align)

    print "rendering..."
    render.run(output_fname, format=format)
    


