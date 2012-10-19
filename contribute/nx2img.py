import networkx as nx
import sys
import os
from imggenerator import SchrodImgGenerator, VirtualX
EDGE_FONTSIZE = 40
def guess_format(fname):
    format = None
    root, ext = os.path.splitext(fname)
    if ext:
        format = ext[1:].lower() 
    
    return format

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

    parser.add_option("-t", "--text",
                      action="store_true", dest="text", default=False,
                      help="Dispaly text instead of 2D image for each node.")

    parser.add_option("-m", "--mcs",
                      action="store_true", dest="mcs", default=False,
                      help="Dispaly maximum common substructure as 2D image for each edge.")

    parser.add_option("-f", "--format",
                      action="store", dest="format", default='svg',
                      choices=["bmp", "jpg", "jpeg", "png", "ppm", "tiff", "xbm", "xpm", "svg"],
                      help="Intermediate image format.")
    
    options, arg = parser.parse_args()
    if len(arg) != 2:
        parser.print_help()
        sys.exit(1)

    input_pickle = arg[0]
    output_fname = arg[1]
    basename, suffix = os.path.splitext(os.path.basename(output_fname))
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
    if not options.text:
        print "generating images for molecules..."
    img_list = []
    curr_dir = os.path.abspath(os.path.curdir)
    for node in D.get_node_list():
        name = node.get_name()
        name = name.strip('"')
        smiles = node.get('SMILES')
        if options.text:
            pass
        else:
            img_fname = '%s/%s_%s.%s'%(curr_dir,basename, name,options.format)
            img_generator.generate(smiles, img_fname)
            img_list.append(img_fname)
            #node.set('image', img_fname)    
            title = node.get('title')
            # use space as label
            #node.set('label', title)
            node.set_label('<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0"><TR><TD><IMG SRC="%s"/></TD></TR><TR><TD>%s</TD></TR></TABLE>>'%(img_fname, title))
            node.set_fontsize(EDGE_FONTSIZE)
    
    simi  = [float(edge.get('similarity')) for edge in D.get_edges()]
    scale = 1.0 / max( simi )    
    for i, edge in enumerate(D.get_edge_list()):
        
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
        
        if options.mcs:
            try:
                smiles = edge.get('SMILES')
            except KeyError, e:
                smiles = None
            if smiles:
                img_fname = "%s/%s_edge_%06d.%s"%(curr_dir,basename, i,options.format)
                img_generator.generate(smiles, img_fname)
                img_list.append(img_fname)
                source_id = edge.get_source()
                dest_id = edge.get_destination()
                source = D.get_node(source_id)[0]
                dest = D.get_node(dest_id)[0]
                source_title = source.get('title')
                dest_title = dest.get('title')            
                edge.set_label('<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0"><TR><TD>%s&nbsp;&lt;--&gt;&nbsp;%s</TD></TR><TR><TD><IMG SRC="%s"/></TD></TR></TABLE>>'%(source_title, dest_title, img_fname))
                edge.set_fontsize(EDGE_FONTSIZE)

    print "rendering..."
    if not format:
        print "unrecognized output format"
        sys.exit(1)
    else:
        print "output format: %s"%format
    
    D.write(output_fname, format=format)
    if not options.save:
        for image in img_list:
            os.remove(image)

