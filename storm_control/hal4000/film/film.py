#!/usr/bin/env python
"""
Handles filming.

This module is responsible for everything related to filming,
including starting and stopping the cameras, saving the frames,
etc..

Much of the logic in the main part of Python2/PyQt4 HAL is 
now located in this module.

Hazen 01/17
"""

import os

from PyQt5 import QtWidgets

import storm_control.sc_library.parameters as params

import storm_control.hal4000.halLib.halMessage as halMessage
import storm_control.hal4000.halLib.halModule as halModule
import storm_control.hal4000.qtdesigner.film_ui as filmUi


class FilmBox(QtWidgets.QGroupBox):
    """
    The UI.
    """
    def __init__(self, parameters = None, **kwds):
        super().__init__(**kwds)
        self.parameters = parameters

        # Add default film parameters.        
        self.parameters.add("acq_mode", params.ParameterSetString("Acquisition mode",
                                                                  "acq_mode",
                                                                  "fixed_length",
                                                                  ["run_till_abort", "fixed_length"]))
        
        self.parameters.add("auto_increment", params.ParameterSetBoolean("Automatically increment movie counter between movies",
                                                                         "auto_increment",
                                                                         True))
        
        self.parameters.add("auto_shutters", params.ParameterSetBoolean("Run shutters during the movie",
                                                                        "auto_shutters",
                                                                        True))
        
        self.parameters.add("filename", params.ParameterString("Current movie file name",
                                                               "filename",
                                                               "movie"))
        
        self.parameters.add("filetype", params.ParameterSetString("Movie file type",
                                                                  "filetype",
                                                                  ".dax",
                                                                  [".dax", ".tif"]))
        
        self.parameters.add("frames", params.ParameterRangeInt("Movie length in frames",
                                                               "frames",
                                                               10,
                                                               1,
                                                               1000000000))
        
        self.parameters.add("want_bell", params.ParameterSetBoolean("Sound bell at the end of long movies",
                                                                    "want_bell",
                                                                    True))
        
        self.parameters.add("dax_big_endian", params.ParameterSetBoolean("Save .dax movies using a big endian format",
                                                                         "dax_big_endian",
                                                                         False))

        # Initial UI configuration.
        self.ui = filmUi.Ui_GroupBox()
        self.ui.setupUi(self)
        
        for extname in self.parameters.getp("extension").getAllowed():
            self.ui.extensionComboBox.addItem(extname)
        
        for typename in self.parameters.getp("filetype").getAllowed():
            self.ui.filetypeComboBox.addItem(typename)

        self.ui.framesText.setText("")
        self.ui.sizeText.setText("")

        self.setDirectory(self.parameters.get("directory"))
        self.setShutters("NA")
        self.newParameters(self.parameters)
        self.updateFilenameLabel()

        # Connect signals
        self.ui.autoIncCheckBox.stateChanged.connect(self.handleAutoInc)
        self.ui.autoShuttersCheckBox.stateChanged.connect(self.handleAutoShutters)
        self.ui.extensionComboBox.currentIndexChanged.connect(self.handleExtension)
        self.ui.filenameEdit.textChanged.connect(self.handleFilename)
        self.ui.filetypeComboBox.currentIndexChanged.connect(self.handleFiletype)
        self.ui.indexSpinBox.valueChanged.connect(self.handleIndex)
        self.ui.lengthSpinBox.valueChanged.connect(self.handleLength)
        self.ui.modeComboBox.currentIndexChanged.connect(self.handleMode)

    def handleAutoInc(self, state):
        self.parameters.set("auto_increment", state)

    def handleAutoShutters(self, state):
        self.parameters.set("auto_shutters", state) 

    def handleExtension(self, index):
        self.parameters.set("extension", self.ui.extensionComboBox.currentText())
        self.updateFilenameLabel()

    def handleFilename(self):
        self.parameters.set("filename", str(self.ui.filenameEdit.displayText()))
        self.updateFilenameLabel()

    def handleFiletype(self, index):
        self.parameters.set("filetype", str(self.ui.filetypeComboBox.currentText()))
        self.updateFilenameLabel()
        
    def handleIndex(self, index):
        self.updateFilenameLabel()

    def handleLength(self, index):
        self.parameters.set("frames", index)

    def handleMode(self, index):
        print("hm", index)
        if (index == 0):
            self.parameters.set("acq_mode", "run_till_abort")
            self.ui.lengthSpinBox.hide()
        else:
            self.parameters.set("acq_mode", "fixed_length")
            self.ui.lengthSpinBox.show()
        
    def newParameters(self, parameters):
        self.ui.autoIncCheckBox.setChecked(parameters.get("auto_increment"))
        self.ui.autoShuttersCheckBox.setChecked(parameters.get("auto_shutters"))
        self.ui.extensionComboBox.setCurrentIndex(self.ui.extensionComboBox.findText(parameters.get("extension")))
        self.ui.filenameEdit.setText(parameters.get("filename"))
        self.ui.filetypeComboBox.setCurrentIndex(self.ui.filetypeComboBox.findText(parameters.get("filetype")))
        self.ui.lengthSpinBox.setValue(parameters.get("frames"))
        
        if (parameters.get("acq_mode") == "run_till_abort"):
            self.ui.modeComboBox.setCurrentIndex(0)
        else:
            self.ui.modeComboBox.setCurrentIndex(1)

    def setDirectory(self, new_directory):
        self.parameters.set("directory", new_directory)
        self.ui.directoryText.setText("  " + new_directory[-30:])

    def setShutters(self, new_shutters):
        self.ui.shuttersText.setText("  " + new_shutters)

    def updateFilenameLabel(self):
        name = self.parameters.get("filename")
        name += "_{0:04d}".format(self.ui.indexSpinBox.value())
        if len(self.parameters.get("extension")) > 0:
            name += "_" + self.parameters.get("extension")
        name += self.parameters.get("filetype")

        self.ui.filenameLabel.setText(name)
        if os.path.exists(os.path.join(self.parameters.get("directory"), name)):
            self.will_overwrite = True
            self.ui.filenameLabel.setStyleSheet("QLabel { color: red}")
        else:
            self.will_overwrite = False
            self.ui.filenameLabel.setStyleSheet("QLabel { color: black}")
        

class Film(halModule.HalModuleBuffered):

    def __init__(self, module_params = None, qt_settings = None, **kwds):
        super().__init__(**kwds)

        self.logfile_fp = open(module_params.get("directory") + "image_log.txt", "a")

        p = module_params.getp("parameters")
        p.add("directory", params.ParameterStringDirectory("Current working directory",
                                                           "directory",
                                                           module_params.get("directory"),
                                                           is_mutable = False,
                                                           is_saved = False))
                
        self.view = FilmBox(parameters = p)

        self.configure_dict = {"ui_order" : 1,
                               "ui_parent" : "hal.containerWidget",
                               "ui_widget" : self.view}

    def cleanUp(self, qt_settings):
        self.logfile_fp.close()
        
    def processMessage(self, message):
        super().processMessage(message)
        if (message.level == 1):
            if (message.m_type == "configure"):
                self.newMessage.emit(halMessage.HalMessage(source = self,
                                                           m_type = "add to ui",
                                                           data = self.configure_dict))


#
# The MIT License
#
# Copyright (c) 2017 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
