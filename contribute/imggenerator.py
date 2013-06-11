import os
import time
import subprocess
import sys
class ImgGenerator:
    """
    A base class to generate 2D structure image file from SMILES
    """
    
    def __init__(self):
        self._supported_formats = {}
    
    def getSupportedFormats(self):
        return self._supported_formats.keys()

        
    def generate(self, smiles, image, format=None):
        """
        Generating 2D image file from SMILES string
        @type smiles: string
        @param smiles: SMILES string representing molecule structure.
        @type image: string
        @param image: The output image file name of 2D structure.
        @type format: string
        @param format: output format of the image file. By default,
        the output format will be deduced from the output image file name.
        """
        if format is None:
            root, ext = os.path.splitext(image)
            if ext is None:
                raise RuntimeError("can not determine image format from suffix")
            format = ext[1:].lower()
        else:
            format = format.lower()
        
        if format not in self._supported_formats.keys():
            raise RuntimeError("format (%s) is not supported"%format)
        
        image_generator = self._supported_formats[format]
        image_generator(smiles, image)
        
    def register(self, format, generator):
        self._supported_formats[format] = generator


class VirtualX(object):
    """
    A class to provide a virtual X server.

    """

    def __init__(self):
        """
        Start up Xvfb and set the DISPLAY environment variable.

        """
        if 'DISPLAY' in os.environ:
            raise RuntimeError("DISPLAY is already defined.")

        self.sleep = time.sleep

        # Start with a high number in the hope that it won't be used. But,
        # check to see if Xvfb startup fails and then keep bumping the
        # display number until it succeeds.

        display_number = 55
        while True:
            self.display = ":%d" % display_number
            self.process = subprocess.Popen(["Xvfb", self.display],
                                            stderr=subprocess.PIPE)
            self.sleep(1)
            if self.process.poll():
                display_number += 1
            else:
                break

        # Make sure that xvfb gets killed if the driver is killed by a
        # SIGTERM.
        import signal
        def kill_xvfb(signum, frame):
            self.terminate()
            print >>sys.stderr, "Terminated"
            # Exit with os._exit; otherwise the unittest runner just
            # catches the SystemExit and continues.
            os._exit(1)
        signal.signal(signal.SIGTERM, kill_xvfb)

        os.environ['DISPLAY'] = self.display

    def terminate(self):
        """
        Kill Xvfb and unset the DISPLAY variable.

        """
        if os and os.environ and 'DISPLAY' in os.environ:
            del os.environ['DISPLAY']
        self.process.terminate()
        # If the process hasn't terminated, sleep for a second and kill it.
        if self.process.returncode is None:
            self.sleep(1)
            self.process.kill()
        self.process = None

    def __del__(self):
        """
        Terminate Xvfb on garbage collection if it wasn't already killed.

        """
        if self.process:
            self.terminate()
try:    
    class SchrodImgGenerator(ImgGenerator):
        """
        A class to generate 2D structure image file from SMILES using Schrodinger toolkit.
        """
        def __init__(self):
        
            ImgGenerator.__init__(self)

            self._load_modules()
        
            self._app = self.QtGui.QApplication([])
            self._model = self.canvas2d.ChmRender2DModel()
            self._renderer = self.canvas2d.Chm2DRenderer(self._model)    

            self.register('bmp', self.smiles2image)
            self.register('jpg', self.smiles2image)
            self.register('jpeg', self.smiles2image)
            self.register('png', self.smiles2image)
            self.register('ppm', self.smiles2image)
            self.register('tiff', self.smiles2image)
            self.register('xbm', self.smiles2image)
            self.register('xpm', self.smiles2image)
            self.register('svg', self.smiles2svg)
        
        def _load_modules(self):
            from PyQt4 import QtGui, QtCore, QtSvg
            from schrodinger.application.canvas.base import ChmLicenseShared_isValid
            from schrodinger.application.canvas.utils import get_license
            from schrodinger.infra import canvas2d
            if not ChmLicenseShared_isValid():
                self.canvas_license = get_license("LICENSE_SHARED")        
            self.QtGui = QtGui
            self.QtCore = QtCore
            self.QtSvg = QtSvg
            self.canvas2d = canvas2d
            virtual_x = None
            if sys.platform != 'win32' and 'DISPLAY' not in os.environ:
                # If there is no DISPLAY variable but there are tests that need X,
                # start Xvfb.
                print "Using VirtualX"
                virtual_x = VirtualX()          
    
        
        def smiles2mol(self, smiles):
            return self.canvas2d.ChmMol.fromSMILES(smiles)    

        def mol2svg(self, mol, svg_fname):
            pic = self._renderer.getQPicture(mol)
            rect = pic.boundingRect()
            svg_gen = self.QtSvg.QSvgGenerator()
            svg_gen.setFileName(svg_fname)
            svg_gen.setSize(self.QtCore.QSize(rect.width(), rect.height()));
            qp = self.QtGui.QPainter(svg_gen) 
            qp.drawPicture(-pic.boundingRect().left(), -pic.boundingRect().top(), pic)
            qp.end()
        
        def smiles2svg(self, smiles, svg_fname):
            mol = self.smiles2mol(smiles)
            self.mol2svg(mol, svg_fname)

    
        def mol2image(self, mol, image_fname): 
            pic = self._renderer.getQPicture(mol)
            rect = pic.boundingRect()
            qimage = self.QtGui.QImage(rect.width(), rect.height(),self.QtGui.QImage.Format_ARGB32)
            qimage.fill(self.QtCore.Qt.white)
            qp = self.QtGui.QPainter(qimage) 
            qp.drawPicture(-pic.boundingRect().left(), -pic.boundingRect().top(), pic)
            qp.end()
            qimage.save(image_fname)
    
        def smiles2image(self, smiles, image_fname):
            mol = self.smiles2mol(smiles)
            self.mol2image(mol, image_fname)

except ImportError, e :                                       
    pass

try:
    import openeye.oechem as oechem
    import openeye.oedepict as oedepict
    class OpeneyeImgGenerator(ImgGenerator):

        """
        A class to generate 2D structure image file from SMILES using Openeye toolkit.
        """
        def __init__(self):
            ImgGenerator.__init__(self)
            self.register('png', self.smiles2svg)
            self.register('svg', self.smiles2svg)
        
        def smiles2mol(self , smiles):

            mol = oechem.OEGraphMol()
            oechem.OEParseSmiles(mol, smiles)
            return mol
        def mol2svg (self, mol , svg_file ):
            oechem.OEAssignImplicitHydrogens(mol)
            oechem.OEAddExplicitHydrogens(mol)
            oedepict.OEPrepareDepiction(mol)             
            width = 1000
            height = 1000
            opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
            pens = opts.GetDefaultBondPen()
            pens.SetLineWidth(10)
            disp = oedepict.OE2DMolDisplay(mol, opts)
            oedepict.OERenderMolecule(svg_file, disp)
        def smiles2svg(self, smiles, svg_file):
            mol = self.smiles2mol(smiles)
            self.mol2svg(mol, svg_file)        
except ImportError, e :                                       
    pass        
