import networkx as nx
import sys
import os
from imggenerator import SchrodImgGenerator, VirtualX

def guess_format(fname):
    format = None
    root, ext = os.path.splitext(fname)
    if ext:
        format = ext[1:].lower() 
    
    return format

def render_smiles(smiles, img_fname):
    img_generator.generate(smiles, img_fname)
    return '<TR><TD><IMG SRC="%s"/></TD></TR>'%img_fname

def render_attr(attr_name, text):
    return render_text("%s: %s"%(attr_name, text))

def render_text(text):
    return "<TR><TD>%s</TD></TR>"%(text)

class DotRender:
    def __init__(self, D, img_generator, 
                 node_attributes, 
                 edge_attributes, 
                 save_image=True,
                 img_format='svg',
                 font_size=40):
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
    
    def _renderNodes(self):
        for i, node in enumerate(self._D.get_node_list()):    
            node_index = i + 1
            node.set_fontsize(self._font_size)   
            name = node.get_name()
            name = name.strip('"')
            render_string = ''
            for attr_name  in self._node_attributes:
                attr = node.get(attr_name)
                if attr:
                    if attr_name == 'SMILES':
                        img_fname = self._getImgFname(name)
                        img_fname = self._getImgFname("node_%06d"%node_index)
                        render_string = render_string + render_smiles(attr, img_fname)
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
                            render_string = render_string + render_smiles(attr, img_fname)
                            self._img_list.append(img_fname)
                        else:
                            render_string = render_string + render_attr(attr_name, attr)
                        
                except KeyError, e:
                    print e
                    pass        
            
            edge.set_label('<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0">%s</TABLE>>'%(render_string)) 
        


    




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
    
    options, arg = parser.parse_args()
    if len(arg) != 2:
        parser.print_help()
        sys.exit(1)

    if len(options.node_attributes) == 0:
        options.node_attributes.append('SMILES')
        options.node_attributes.append('title')
        
    if len(options.edge_attributes) == 0:
        options.node_attributes.append('SMILES')
  
        
    input_pickle = arg[0]
    output_fname = arg[1]
    
    format = guess_format(output_fname)
    if options.format == 'svg' and format == 'svg' and not options.save:
        print "WARNING:the output svg file depends on intermediate svg files, please use -s option."
        sys.exit(0)

    G = nx.read_gpickle(input_pickle)
    D = nx.to_pydot(G)
    
    virtual_x = None
    if sys.platform != 'win32' and 'DISPLAY' not in os.environ:
        # If there is no DISPLAY variable but there are tests that need X,
        # start Xvfb.
        print "Using VirtualX"
        virtual_x = VirtualX()    
    img_generator = SchrodImgGenerator()

    render = DotRender(D, img_generator, 
                       options.node_attributes,
                       options.edge_attributes,
                       save_image=options.save,
                       img_format=options.format)

    print "rendering..."
    render.run(output_fname, format=format)
    


