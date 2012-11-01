from pyquery import PyQuery as pq
from lxml import etree
import sys
from pyquery.pyquery import fromstring
import os
import argparse
import xml.etree.ElementTree as ET
ET.register_namespace('', "http://www.w3.org/2000/svg")
ET.register_namespace('xlink', "http://www.w3.org/1999/xlink")

HREF = '{http://www.w3.org/1999/xlink}href'

# pyquery has to use it to workaround namespace issues in dealing with xml
namespaces = {"s":"http://www.w3.org/2000/svg", 
               "x":"http://www.w3.org/1999/xlink"}

html_template=\
"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>LeadOptMap</title>

    <!-- Le styles -->
    <link href="docs/assets/css/bootstrap.css" rel="stylesheet">
    <link href="docs/assets/css/docs.css" rel="stylesheet">
    <link href="docs/assets/js/google-code-prettify/prettify.css" rel="stylesheet">

    <!-- CSS just for the tests page -->
    <!-- <link href="less/tests/css-tests.css" rel="stylesheet"> -->

    <script type="text/javascript" src="js/jquery-1.8.2.js"></script>
    <script type="text/javascript" src="js/dracula_graph.js"></script>
    <script type="text/javascript" src="js/sprintf-0.6.js"></script>
    <script type="text/javascript" src="js/leadoptmap.js"></script>

</head>
<body>

<!-- Masthead
================================================== -->
<header class="jumbotron subhead" id="overview">
    <div class="container">
        <h1>Lead Optimization Map</h1>
    </div>
</header>


<div class="bs-docs-canvas">

<div class="container">

<!-- Tables
================================================== -->


<div class="row">
    <h2>Mutation Pair and Maximum Common Substructure</h2>
    <table id="mt"  class="table table-bordered">
        <colgroup id="mt_colorgroup">
            <col class="col1">
            <col class="col2">
        </colgroup>
        <thead>
        <tr  id="mt_thead_tr">
            <th>Ligand1</th>
            <th>Ligand2</th>

        </tr>
        </thead>
        <tbody id="mt_tbody">

        </tbody>
        <tfoot>
        <tr id="mt_dropdown">
            <td>
                <select id="mol1" class="span4">Ligand1 </select>
            </td>
            <td>
                <select id="mol2" class="span4">Ligand2</select>
            </td>
        </tr>
        </tfoot>
    </table>
    <div class="row-fluid">
        <button id="addEdgeButton" class="btn btn-success">Add Edge</button>
        <button id="addAllEdgesButton" class="btn btn-success">Add All Edges</button>
        <button id="clearButton" class="btn btn-success">Clear</button>

    </div>
</div><!--/row-->

<div class="row">
    <h2>Control Using Jquery</h2>
    <div class="row-fluid">
        <div  class="span12">
            <textarea id="jquery_command" class="span12"></textarea>
        </div>
    </div>

    <!--do not put them inside form, otherwise, the refresh will break-->
    <div class="row-fluid">
        <button id="runButton" class="btn btn-success">Run</button>
        <button id="showAllButton" class="btn btn-success">Show All</button>
        <button id="hideAllButton"  class="btn btn-info">Hide All</button>
    </div>

</div><!--/row-->
<div id="status" class="row"></div>
<br>

<!-- SVG
================================================== -->

<div id="map"  class="container">
%s
</div>
</div>
</body>
</html>
"""
def parse_svg(filename):
    context = open(filename).read()
    #tree = ET.parse(filename)
    #xml = fromstring(context, parser='xml', custom_parser=ET.parse)
    
    #d = pq([tree.getroot()])   
    d = pq(filename=filename, parser='xml')
    return d


def replace_svg(m, image, attrs):
    href = m(image).attr(HREF) 
    c = parse_svg(href)
    c('s|svg titlei, desc', namespaces=namespaces).remove()
    
    x = image.attrib['x']
    y = image.attrib['y']
    attrs['transform']="translate(%s, %s)"%(x, y)
    attrs['href']='img/'+href
    
    attr_str = ' '.join(['%s="%s"'%(key,val) for key, val in attrs.items()])
    html = c('s|svg', namespaces=namespaces).html()
    html = '<g %s> %s </g>'%(attr_str, html)
    return m(image).replaceWith(html)    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
"""
A script to compare fepio_fep block and associated ffio block between two cms files.
Note that, the coordinates in these two files must match each other.

""")
    parser.add_argument('input', action='store',
                       help='maestro file')
    parser.add_argument('output', action='store',
                       help='maestro file')
    parser.add_argument('-c', '--combine', dest='combine', action='store_true',
                       default=False, help='combine various svg files into html')
    
    parser.add_argument('-t', '--target', dest='target', action='store',
                       default=None, help='')    

    args = parser.parse_args()
            
    input_fname = args.input
    output_fname = args.output 

    
    
    # parse main svg file
    m = parse_svg(input_fname)
    graph = m('s|svg > s|g', namespaces=namespaces)
    transform = graph.attr('transform')
    scale = 'scale(1 1)' 
    if scale in transform:
        transform = 'scale(0.1 0.1)' + transform[len(scale):]    
        graph.attr('transform', transform)
    
    nodes = {}
    edges = []     

    all_images = []
    for node in m('s|g[id^="node"]', namespaces=namespaces):
        title = m(node)('s|title', namespaces=namespaces)[0]
        
        text = m(node)('s|text', namespaces=namespaces)[0]

        nid = title.text.strip()
        name = text.text.strip()
        if name.startswith('title:'):
            name = name[len('title:'):]
        
        title.text = name
        node.attrib['name'] = name
        nodes[nid] = node
        image = m(node)('s|image', namespaces=namespaces)[0]
            
        
        href = m(image).attr(HREF) 
        all_images.append(href)
        if args.combine:
            attrs = {'class':'svg',
                     'node_id':node.attrib['id']
                    }
            replace_svg(m, image, attrs)        
        else:
            
            m(image).attr(HREF, 'img/'+href)
            
        
        
    for edge in m('s|g[id^="edge"]', namespaces=namespaces):
        title = m(edge)('s|title', namespaces=namespaces)[0]
        
        text = m(edge)('s|text', namespaces=namespaces)[0]

        nids = title.text.split('--')
        if len(nids) != 2:
            raise Exception("")
        source_node = nodes[nids[0]]
        target_node = nodes[nids[1]]
        edge.attrib['source'] = source_node.attrib['id']
        edge.attrib['target'] = target_node.attrib['id']
        source_name = source_node.attrib['name']
        target_name = target_node.attrib['name']        
        
        title.text='%s:%s'%(source_name, target_name)

        image_list = m(edge)('s|image', namespaces=namespaces)
        for image in image_list:
            
            href = m(image).attr(HREF)
            all_images.append(href)
            prefix, suffix = os.path.splitext(href)
            basename = os.path.basename(prefix)
            tokens = basename.split('_')
            if len(tokens) > 2:
                mcs_type = tokens[-1]
                node_id = tokens[-2]
                if not node_id.startswith('node'):
                    raise Exception("invalid node_id:%s"%node_id)
            
            if args.combine:
                attrs = {'class':'svg',
                         'node_id':node_id
                        }                
                replace_svg(m, image, attrs)
            else:
                image.attrib[HREF] = 'img/'+href
                image.attrib['mcs_type'] = mcs_type
                image.attrib['node_id'] =node_id
                
    
    if args.target:
        output_fname = os.path.join(args.target, os.path.basename(output_fname))
    if output_fname.endswith('.htm') or output_fname.endswith('.html'):   
        open(output_fname, 'w').write(html_template%str(m))
    else:
        open(output_fname, 'w').write(str(m))
    
    if args.target:
        import shutil
        target = os.path.join(args.target, 'img')
        for img in all_images:
            shutil.copy(img, target)
