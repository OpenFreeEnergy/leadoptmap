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

if __name__ == '__main__':
    usage = \
    """usage: %prog [options] <input.pkl> <output.png>
    Generating an image file from a networkx graph pickle file. 
    Each node must have a SMILES attribute that represents the associated molecule.
    """
    from optparse import OptionParser
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--save",
                      action="store_true", dest="save", default=False,
                      help="Save intermediate images for individual molecule.")

    parser.add_option("-t", "--text",
                      action="store_true", dest="text", default=False,
                      help="Dispaly text instead of 2D image.")

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
            img_fname = '%s/%s.svg'%(curr_dir,name)
            img_generator.generate(smiles, img_fname)
            img_list.append(img_fname)
            node.set('image', img_fname)        
            # use space as label
            node.set('label', ' ')
        

        
    print "rending..."
    format = guess_format(output_fname)
    if not format:
        print "unrecognized output format"
        sys.exit(1)
    else:
        print "output format: %s"%format
    
    D.write(output_fname, format=format)
    if not options.save:
        for image in img_list:
            os.remove(image)

