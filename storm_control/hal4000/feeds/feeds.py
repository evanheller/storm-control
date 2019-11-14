#!/usr/bin/env python
"""
This module enables the processing of camera frame(s) with
operations like averaging, slicing, etc..

It is also responsible for keeping tracking of how many
different cameras / feeds are available for each parameter
file, whether the cameras / feeds should be saved when
filming and what extension to use when saving.

Hazen 03/17
"""

import copy
import numpy

from PyQt5 import QtCore, QtWidgets

import storm_control.sc_library.hdebug as hdebug
import storm_control.sc_library.halExceptions as halExceptions
import storm_control.sc_library.parameters as params

import storm_control.hal4000.camera.frame as frame
import storm_control.hal4000.camera.cameraFunctionality as cameraFunctionality
import storm_control.hal4000.halLib.halMessage as halMessage
import storm_control.hal4000.halLib.halMessageBox as hmb
import storm_control.hal4000.halLib.halModule as halModule

# Added by EH for GUI
import storm_control.hal4000.halLib.halDialog as halDialog
import storm_control.hal4000.qtdesigner.timelapse_ui as timelapseUi
import storm_control.hal4000.illumination.xmlParser as xmlParser


def checkParameters(parameters):
    """
    Checks parameters to verify that there won't be any errors
    when we actually try to create the feeds.

    Throw an exception if there are any problems.
    """
    if not parameters.has("feeds"):
        return

    # Check the feed parameters. For now all we are doing is verifying that
    # the feed ROI area is a multiple of 4.
    #
    feed_parameters = parameters.get("feeds")
    for feed_name in feed_parameters.getAttrs():
        fp = feed_parameters.get(feed_name)
        cp = parameters.get(fp.get("source"))

        x_start = fp.get("x_start", 1)
        x_end = fp.get("x_end", cp.get("x_pixels"))
        x_pixels = x_end - x_start + 1

        y_start = fp.get("y_start", 1)
        y_end = fp.get("y_end", cp.get("y_pixels"))
        y_pixels = y_end - y_start + 1
        
        # Check that the feed size is a multiple of 4 in x.
        if not ((x_pixels % 4) == 0):
            raise FeedException("The x size of the feed ROI must be a multiple of 4 in " + feed_name)


class FeedException(halExceptions.HalException):
    pass


class FeedFunctionality(cameraFunctionality.CameraFunctionality):
    """
    Feed functionality in a form that other modules can interact with. These have
    a camera functionality which they are interacting with to create the feed.

    For the most part this just passes through information from the underlying
    camera functionality.

    Some functionality is explicitly blocked so we get an error if we accidentally
    try and use this exactly like a camera functionality.
    """
    def __init__(self, feed_name = None, **kwds):
        super().__init__(**kwds)
        self.cam_fn = None
        self.feed_name = feed_name
        self.feed_parameters = self.parameters
        self.frame_number = 0
        self.frame_slice = None
        self.number_connections = 0
        self.x_pixels = 0
        self.y_pixels = 0

        # We're going to change some of the parameters in this class, but we don't
        # want the new values to appear in the editor or when saved. In order to
        # do this we change self.parameters to be a copy of self.feed_parameters.
        #
        self.parameters = self.feed_parameters.copy()

    def connectCameraFunctionality(self):
        """
        Connect the feed to it's camera functionality
        """
        # sanity check.
        assert(self.number_connections == 0)
        self.number_connections += 1
        
        self.cam_fn.newFrame.connect(self.handleNewFrame)
        self.cam_fn.started.connect(self.handleStarted)
        self.cam_fn.stopped.connect(self.handleStopped)

    def disconnectCameraFunctionality(self):
        """
        Disconnect the feed from it's camera functionality.
        """
        # sanity check.
        assert(self.number_connections == 1)
        self.number_connections += 1
        
        if self.cam_fn is not None:
            self.cam_fn.newFrame.disconnect(self.handleNewFrame)
            self.cam_fn.started.disconnect(self.handleStarted)
            self.cam_fn.stopped.disconnect(self.handleStopped)

    def getCameraFunctionality(self):
        """
        Return the camera functionality this feed is using.
        """
        return self.cam_fn

    def getFeedName(self):
        """
        Return the name of the feed (as specified in the XML file).
        """
        return self.feed_name

    def handleNewFrame(self, new_frame):
        sliced_data = self.sliceFrame(new_frame)
        self.newFrame.emit(frame.Frame(sliced_data,
                                       new_frame.frame_number,
                                       self.x_pixels,
                                       self.y_pixels,
                                       self.camera_name))

    def updateParameters(self, parameters):
        params.copyParametersReplace("", self.parameters, parameters)
        #self.parameters.copyParameters(parameters)
            
    def handleStarted(self):
        self.started.emit()

    def handleStopped(self):
        self.stopped.emit()

    def hasEMCCD(self):
        assert False
        
    def hasPreamp(self):
        assert False

    def hasShutter(self):
        assert False

    def hasTemperature(self):
        assert False

    def haveCameraFunctionality(self):
        return (not self.cam_fn is None)
        
    def isCamera(self):
        return False

    def isMaster(self):
        return False

    def reset(self):
        self.frame_number = 0

    def setCameraFunctionality(self, camera_functionality):
        self.cam_fn = camera_functionality

        # We are changing some of the parameters here so that the feed will
        # be displayed properly.

        # The assumption here is that x_start, x_end and x_pixels are all in
        # units of binned pixels. This is also what we assume with the camera.
        #
        # Also, the initial values for x_start and x_end will be 1, if they
        # were not specified in the parameters file.
        #
        p = self.parameters

        # Figure out if we need to slice.
        if (p.get("x_end") == 1):
            p.setv("x_end", self.cam_fn.getParameter("x_end"))
        if (p.get("y_end") == 1):
            p.setv("y_end", self.cam_fn.getParameter("y_end"))

        p.set("x_pixels", p.get("x_end") - p.get("x_start") + 1)
        p.set("y_pixels", p.get("y_end") - p.get("y_start") + 1)
        p.set("bytes_per_frame", 2 * p.get("x_pixels") * p.get("y_pixels"))

        self.x_pixels = p.get("x_pixels")
        self.y_pixels = p.get("y_pixels")

        if (p.get("x_pixels") != self.cam_fn.getParameter("x_pixels")) or\
           (p.get("y_pixels") != self.cam_fn.getParameter("y_pixels")):
            self.frame_slice  = (slice(p.get("y_start") - 1, p.get("y_end")),
                                 slice(p.get("x_start") - 1, p.get("x_end")))

        # Adjust / add parameters from the camera so that the feed will be
        # displayed properly.
        p.add(self.cam_fn.parameters.getp("x_chip").copy())
        p.add(self.cam_fn.parameters.getp("y_chip").copy())

        p.setv("x_start", self.cam_fn.getParameter("x_start") + p.get("x_start") - 1)
        p.setv("x_end", p.get("x_start") + p.get("x_pixels") - 1)
        
        p.setv("y_start", self.cam_fn.getParameter("y_start") + p.get("y_start") - 1)
        p.setv("y_end", p.get("y_start") + p.get("y_pixels") - 1)

        # Set the maximums of the original versions of the parameters. These
        # are what will be passed to the parameters editor.
        #
        self.feed_parameters.getp("x_start").setMaximum(self.cam_fn.getParameter("x_pixels"))
        self.feed_parameters.getp("x_end").setMaximum(self.cam_fn.getParameter("x_pixels"))
        self.feed_parameters.getp("y_start").setMaximum(self.cam_fn.getParameter("y_pixels"))
        self.feed_parameters.getp("y_end").setMaximum(self.cam_fn.getParameter("y_pixels"))

        # Add some of other parameters that we need to behave like a camera functionality. These
        # are just duplicates from the corresponding camera.
        for pname in ["default_max", "default_min", "flip_horizontal", "flip_vertical",
                      "fps", "max_intensity", "transpose", "x_bin", "y_bin"]:
            p.add(self.cam_fn.parameters.getp(pname).copy())

        # Connect camera functionality signals. We just pass most of
        # these through.
        self.connectCameraFunctionality()

    def sliceFrame(self, new_frame):
        """
        Slices out a part of the frame based on self.frame_slice.
        """
        if self.frame_slice is None:
            return new_frame.np_data
        else:
            w = new_frame.image_x
            h = new_frame.image_y
            sliced_frame = numpy.reshape(new_frame.np_data, (h,w))[self.frame_slice]
            return numpy.ascontiguousarray(sliced_frame)

    def toggleShutter(self):
        assert False


class FeedFunctionalityAverage(FeedFunctionality):
    """
    The feed functionality for averaging frames together.
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)

        self.average_frame = None
        self.counts = 0
        self.frames_to_average = self.parameters.get("frames_to_average")

    def handleNewFrame(self, new_frame):
        sliced_data = self.sliceFrame(new_frame)

        if self.average_frame is None:
            self.average_frame = sliced_data.astype(numpy.uint32)
        else:
            self.average_frame += sliced_data
        self.counts += 1

        if (self.counts == self.frames_to_average):
            average_frame = self.average_frame/self.frames_to_average
            self.newFrame.emit(frame.Frame(average_frame.astype(numpy.uint16),
                                           self.frame_number,
                                           self.x_pixels,
                                           self.y_pixels,
                                           self.camera_name))
            self.average_frame = None
            self.counts = 0
            self.frame_number += 1

    def reset(self):
        super().reset()
        self.average_frame = None
        self.counts = 0
        
    
class FeedFunctionalityInterval(FeedFunctionality):
    """
    The feed functionality for picking out a sub-set of the frames.
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)

        temp = self.parameters.get("capture_frames")
        self.capture_frames = list(map(int, temp.split(",")))
        self.cycle_length = self.parameters.get("cycle_length")

    def handleNewFrame(self, new_frame):
        sliced_data = self.sliceFrame(new_frame)
        
        if (new_frame.frame_number % self.cycle_length) in self.capture_frames:
            self.newFrame.emit(frame.Frame(sliced_data,
                                           self.frame_number,
                                           self.x_pixels,
                                           self.y_pixels,
                                           self.camera_name))
            self.frame_number += 1

    def updateParameters(self, parameters):
        super().updateParameters(parameters)
        self.cycle_length = self.parameters.get("cycle_length")


class FeedFunctionalitySlice(FeedFunctionality):
    """
    The feed functionality for slicing out sub-sets of frames.
    """
    pass

        
class FeedController(object):
    """
    Feed controller.
    """
    def __init__(self, parameters = None, **kwds):
        """
        parameters - This is just the 'feed' section of the parameters.
        """
        super().__init__(**kwds)

        self.feeds = {}
        if parameters is None:
            return

        # Create the feeds.
        self.parameters = parameters
        for feed_name in self.parameters.getAttrs():
            file_params = self.parameters.get(feed_name)
            
            # Create default feed parameters.
            max_value = 100000
            feed_params = params.StormXMLObject()

            # Feeds are saved with their name as the extension.
            feed_params.add(params.ParameterString(name = "extension",
                                                   value = feed_name,
                                                   is_mutable = True))
            
            feed_params.add(params.ParameterString(name = "feed_type",
                                                   value = "",
                                                   is_mutable = False))

            feed_params.add(params.ParameterSetBoolean(name = "saved",
                                                       value = False))

            # This is the camera that drives the feed.
            feed_params.add(params.ParameterString(name = "source",
                                                   value = "",
                                                   is_mutable = False))
            
            feed_params.add(params.ParameterRangeInt(description = "AOI X start.",
                                                     name = "x_start",
                                                     value = 1,
                                                     min_value = 1,
                                                     max_value = max_value))

            feed_params.add(params.ParameterRangeInt(description = "AOI X end.",
                                                     name = "x_end",
                                                     value = 1,
                                                     min_value = 1,
                                                     max_value = max_value))

            feed_params.add(params.ParameterRangeInt(description = "AOI Y start.",
                                                     name = "y_start",
                                                     value = 1,
                                                     min_value = 1,
                                                     max_value = max_value))
            
            feed_params.add(params.ParameterRangeInt(description = "AOI Y end.",
                                                     name = "y_end",
                                                     value = 1,
                                                     min_value = 1,
                                                     max_value = max_value))

            # Figure out what type of feed this is.
            fclass = None
            feed_type = file_params.get("feed_type")
            
            if (feed_type == "average"):
                fclass = FeedFunctionalityAverage
                
                feed_params.add(params.ParameterInt(description = "Number of frames to average.",
                                                    name = "frames_to_average",
                                                    value = 1))
                            
            elif (feed_type == "interval"):
                fclass = FeedFunctionalityInterval

                feed_params.add(params.ParameterInt(description = "Interval cycle length.",
                                                    name = "cycle_length",
                                                    value = 1))
                
                feed_params.add(params.ParameterCustom(description = "Frames to capture.",
                                                       name = "capture_frames",
                                                       value = "1"))

            elif (feed_type == "slice"):
                fclass = FeedFunctionalitySlice
            else:
                raise FeedException("Unknown feed type '" + feed_type + "' in feed '" + feed_name + "'")

            # Update with values from the parameters file. Depending on the parameters
            # file it might include parameters that we don't have and which we silently
            # ignore.
            #
            for attr in file_params.getAttrs():
                if feed_params.has(attr):
                    feed_params.setv(attr, file_params.get(attr))

            # Replace the values in the parameters that were read from a file with these values.
            self.parameters.addSubSection(feed_name, feed_params, overwrite = True)

            camera_name = feed_params.get("source") + "." + feed_name
            #hmb.halMessageBoxInfo(camera_name)
            
            self.feeds[camera_name] = fclass(feed_name = feed_name,
                                             camera_name = camera_name,
                                             parameters = feed_params)

    def allFeedsFunctional(self):
        for feed in self.getFeeds():
            if not feed.haveCameraFunctionality():
                return False
        return True

    def disconnectFeeds(self):
        """
        Disconnect the feeds from their camera functionalities.
        """
        for feed in self.getFeeds():
            feed.disconnectCameraFunctionality()

    def getFeed(self, feed_name):
        return self.feeds[feed_name]
        
    def getFeedNames(self):
        return list(self.feeds.keys())

    def getFeeds(self):
        return list(self.feeds.values())
    
    def getParameters(self):
        return self.parameters

    def updateParameters(self, parameters):
        params.copyParametersReplace("", self.parameters, parameters)
        
        # Propagate changes to the individual feeds
        for feed in self.getFeeds():
            feed.updateParameters(parameters.get(feed.getFeedName()))

          
    def resetFeeds(self):
        for feed in self.getFeeds():
            feed.reset()
            

class FeedsView(halDialog.HalDialog):
    editParameters = QtCore.pyqtSignal()
    
    def __init__(self, configuration = None, **kwds):
        """
        This initializes things and sets up the UI 
        of the power control dialog box.
        """
        super().__init__(**kwds)

        self.channels = []
        self.configuration = configuration
        self.directory = None
        self.tparams = params.StormXMLObject()
        self.parameters = params.StormXMLObject()
        self.timing_functionality = None
        self.which_checked = []
        self.ilm_functionality = None
        self.exposuretime = None
        self.interval = None

        # UI setup
        self.ui = timelapseUi.Ui_Dialog()
        self.ui.setupUi(self)

        #self.ui.textInterval. 
        self.ui.comboIntervalType.currentTextChanged.connect(self.handleIntervalTypeChange)
        self.ui.sliderInterval.valueChanged.connect(self.handleSlider)
        self.ui.textInterval.textChanged.connect(self.handleIntervalChange)
        #self.ui.textInterval.editingFinished.connect(self.handleIntervalChange)

        # Some initial values
        self.tparams.add(params.ParameterString(name = "extension",
                                               value = "interval",
                                               is_mutable = True))
        
        self.tparams.add(params.ParameterString(name = "feed_type",
                                               value = "interval",
                                               is_mutable = False))

        self.tparams.add(params.ParameterSetBoolean(name = "saved",
                                                   value = False))

        # This is the camera that drives the feed.
        self.tparams.add(params.ParameterString(name = "source",
                                               value = "camera1",
                                               is_mutable = False))

        self.tparams.add(params.ParameterInt(description = "Interval cycle length.",
                                                    name = "cycle_length",
                                                    value = 1))
                
        self.tparams.add(params.ParameterCustom(description = "Frames to capture.",
                                                       name = "capture_frames",
                                                       value = "1"))


        #self.ui.progressionsCheckBox.stateChanged.connect(self.handleProgressionsCheck)

    def setIlmFunctionality(self, functionality):
        """
        Configure the illumination channels using the illumination functionality.
        """
        self.ilm_functionality = functionality
        channels = self.ilm_functionality.getChannelNames()

        parent = self.ui.gridWidget
        layout = parent.layout()
        for i, channel in enumerate(channels):
            self.which_checked.append(False)


            # check box with available laser channels
            channel_active_check_box = QtWidgets.QCheckBox(parent)
            channel_active_check_box.setText(channel)
            layout.addWidget(channel_active_check_box, 1, i+1)

            self.channels.append(channel_active_check_box)

    def setCamFunctionality(self, functionality):
        self.exposuretime = functionality.getParameter("exposure_time")


    def getParameters(self):
        return self.parameters

    def handleSlider(self, value):
        self.ui.textInterval.setText(str(value))

    def handleIntervalChange(self):
        # Handle keeping the interval in units of the current exposure
        cexp = self.exposuretime * 1000
        cval = float(self.ui.textInterval.text())
        fval = float
        
        if self.ui.comboIntervalType.currentText() == 'ms':
            fval = cval - (cval % cexp)
            self.interval = fval
            self.ui.textInterval.setText(str(fval))

        elif self.ui.comboIntervalType.currentText() == 's':
            cval = cval*1000
            fval = cval - (cval % cexp) 
            self.interval = fval
            self.ui.textInterval.setText(str(fval/1000))

        elif self.ui.comboIntervalType.currentText() == 'min':
            cval = cval*1000*60
            fval = cval - (cval % cexp) 
            self.interval = fval
            self.ui.textInterval.setText(str(fval/1000/60))

        elif self.ui.comboIntervalType.currentText() == 'h':
            cval = cval*1000*60*60 
            fval = cval - (cval % cexp) 
            self.interval = fval
            self.ui.textInterval.setText(str(fval/1000/60/60))


        # See if we can't create some parameters out of this
        self.tparams.set("cycle_length", self.interval / cexp+1 ) 
        self.parameters.addSubSection("interval", self.tparams, overwrite = True)
        self.editParameters.emit() 



    def handleIntervalTypeChange(self):
        if self.ui.comboIntervalType.currentText() == 'ms':
            self.ui.sliderInterval.setMaximum(999)

        elif self.ui.comboIntervalType.currentText() == 's':
            self.ui.sliderInterval.setMaximum(600)
            
        elif self.ui.comboIntervalType.currentText() == 'min':
            self.ui.sliderInterval.setMaximum(60)

        else:
            self.ui.sliderInterval.setMaximum(24)

        self.handleIntervalChange()
            

    def startFilm(self):
        
        for i, channel in enumerate(self.channels):
            self.which_checked[i] = False

            if channel.isChecked():
                self.which_checked[i] = True
                

        return [self.which_checked, self.tparams.get("cycle_length")]
                 


class Feeds(halModule.HalModule):
    """
    Feeds controller.
    """
    def __init__(self, module_params = None, qt_settings = None, **kwds):
        super().__init__(**kwds)
        self.camera_names = []
        self.feed_controller = None
        self.feed_names = []
        self.full_params = None
        
        # This message comes from the display.display when it creates a new
        # viewer.
        halMessage.addMessage("get feed names",
                              validator = {"data" : {"extra data" : [False, str]},
                                           "resp" : {"feed names" : [True, list]}})

        halMessage.addMessage("feedparams",
                validator = {"data" : {"parameters" : [True, params.StormXMLObject]},
                                 "resp" : None})

        halMessage.addMessage("new shutters",
                validator = {"data" : {"shutterinfo" : [True, xmlParser.ShuttersInfo],
                                        "waveforms": [True, list],
                                        "oversampling": [True, int]},
                                 "resp" : None})


        configuration = module_params.get("configuration")
        self.ilm_fn_name = configuration.get("illumination_functionality")

        self.view = FeedsView(configuration = configuration,
                                     module_name = self.module_name)
        self.view.halDialogInit(qt_settings,
                                module_params.get("setup_name") + " feeds control")

        self.view.editParameters.connect(self.handleEditParameters)
        

        
    def handleEditParameters(self):
        params = self.view.getParameters()
        checkParameters(params)

        # Case of brand new feed
        if params is not None and self.feed_controller is None:
            self.sendMessage(halMessage.HalMessage(m_type = "initial parameters",
                                                      data = {"parameters" : params }))

            self.feed_controller = FeedController(parameters = params)
            
            # Add camera functionality to all of them
            for feed in self.feed_controller.getFeeds():
                self.feed_names.append(feed.getCameraName())
                self.sendMessage(halMessage.HalMessage(m_type = "get functionality",
                                                       data = {"name" : feed.getParameter("source"),
                                                               "extra data" : feed.getCameraName()}))
       
        else:
            # Pass updated settings onto Settings module
            self.feed_controller.updateParameters(params)

            for feed in self.feed_controller.getFeeds():
                self.sendMessage(halMessage.HalMessage(m_type = "feedparams",
                    data = {"parameters" : params.copy()}))


          
    def broadcastCurrentFeeds(self):
        """
        Send a 'configuration' message with the current feed names.

        film.film uses this message to know what all the feeds are.

        display.cameraDisplay uses this message to populate the feed chooser
        combobox.
        """
        props = {"feed names" : self.feed_names}
        self.sendMessage(halMessage.HalMessage(m_type = "configuration",
                                               data = {"properties" : props}))

    def handleResponse(self, message, response):

        if self.feed_controller is not None:
            if message.isType("get functionality"):
                feed = self.feed_controller.getFeed(message.getData()["extra data"])
                feed.setCameraFunctionality(response.getData()["functionality"])
                self.broadcastCurrentFeeds()

                              
            # If we have camera functionality for all the feeds then it is safe to
            # broadcast the new feed information.
            #
            #if self.feed_controller.allFeedsFunctional():
            #    self.broadcastCurrentFeeds()

                # And we are done with the parameter change. Send the current
                # state of the parameters.
            #    self.sendMessage(halMessage.HalMessage(m_type = "parameters changed",
            #                                          data = {"new parameters" : self.feed_controller.getParameters().copy()}))

        else:
            if message.getData()["name"] == "camera1" :
                self.view.setCamFunctionality(response.getData()["functionality"])
                  # Try saving the cam functionality for later
                self.cam_fun = response.getData()["functionality"]

            else: 
                self.view.setIlmFunctionality(response.getData()["functionality"])

    def processMessage(self, message):

        if message.isType("configuration"):
            # Look for message about the cameras in this setup.
            if ("is camera" in message.getData()["properties"]):
                m_data = message.getData()["properties"]
                if m_data["is camera"]:
                    self.camera_names.append(m_data["module name"])

        elif message.isType("configure1"):
            # Let the settings.settings module know that it needs
            # to wait for us during a parameter change.
            self.sendMessage(halMessage.HalMessage(m_type = "wait for",
                                                  data = {"module names" : ["settings"]}))


            ## 
            # To set up a feed editor dialog without a feeds param already existing,
            # we want to get some information about the illumination available
            # and current camera settings... Adding another message for this so we 
            # don't perturb the existing "get functionality" that was meant to 
            # existing feeds

            # Get information about current illumination settings
            self.sendMessage(halMessage.HalMessage(m_type = "get functionality",
                                                   data = {"name" : self.ilm_fn_name}))

            # Get info about current camera
            # Hardcoded for the moment
            self.sendMessage(halMessage.HalMessage(m_type = "get functionality",
                                                   data = {"name" : "camera1"}))


            self.sendMessage(halMessage.HalMessage(m_type = "add to menu",
                                                   data = {"item name" : "Timelapse",
                                                           "item data" : "timelapse"}))
                                                             
          
        elif message.isType("configure2"):
            self.feed_names = copy.copy(self.camera_names)
            self.broadcastCurrentFeeds()

        # Send film settings along so we can generate a shutter sequence
        elif message.isType("start film"):
            [channels, frames] = self.view.startFilm()
            
            # Attempt building shutter sequence on-the-fly
            oversampling = 100
            waveforms = []
            on =  int( float(0) * float(oversampling))
            off = int( float(1) * float(oversampling))

            color_data = []
            for i in range(frames):
                color_data.append(None)

            for i, c in enumerate(channels):
                waveforms.append(numpy.zeros(frames * oversampling))
                if c is True:
                    waveform = waveforms[i]
                    j = on
                    while j < off:
                        waveform[j] = 1
                        j += 1

            # Propagate the shutters info to other modules
            self.sendMessage(halMessage.HalMessage(m_type = "new shutters",
                data = {"shutterinfo" : xmlParser.ShuttersInfo(color_data= color_data, frames = frames), "waveforms": waveforms, "oversampling": oversampling}))


        elif message.isType("get functionality"):
            if self.feed_controller is not None:
                feed_name = message.getData()["name"]
                if feed_name in self.feed_controller.getFeedNames():
                    message.addResponse(halMessage.HalMessageResponse(source = self.module_name,
                                                                      data = {"functionality" : self.feed_controller.getFeed(feed_name)}))

        elif message.isType("get feed names"):
            message.addResponse(halMessage.HalMessageResponse(source = self.module_name,
                                                              data = {"feed names" : self.feed_names}))
            
        elif message.isType("new parameters"):
            params = message.getData()["parameters"]

            # If WE send the new params method, we don't expect the full list to be here
            #checkParameters(params)
            if self.feed_controller is not None:
                self.feed_controller.disconnectFeeds()
                message.addResponse(halMessage.HalMessageResponse(source = self.module_name,
                                                                  data = {"old parameters" : self.feed_controller.getParameters().copy()}))
                self.feed_controller = None

            if params.has("feeds"):
                self.feed_controller = FeedController(parameters = params.get("feeds"))
            
        elif message.isType("updated parameters"):
            self.feed_names = copy.copy(self.camera_names)

            if self.feed_controller is not None:
                for feed in self.feed_controller.getFeeds():
                    self.feed_names.append(feed.getCameraName())
                    self.sendMessage(halMessage.HalMessage(m_type = "get functionality",
                                                           data = {"name" : feed.getParameter("source"),
                                                                   "extra data" : feed.getCameraName()}))

                # Try moving this here from get_functionality to handle switching between parameters
                if self.feed_controller.allFeedsFunctional():
                    self.broadcastCurrentFeeds()

                # And we are done with the parameter change. Send the current
                # state of the parameters.
                self.sendMessage(halMessage.HalMessage(m_type = "parameters changed",
                                                      data = {"new parameters" : self.feed_controller.getParameters().copy()}))

            else:
                self.broadcastCurrentFeeds()
                self.sendMessage(halMessage.HalMessage(m_type = "parameters changed"))

        elif message.isType("start film"):
            if self.feed_controller is not None:
                self.feed_controller.resetFeeds()

        # Added option to show
        elif message.isType("show"):
            if (message.getData()["show"] == "timelapse"):
                self.view.show()
        
        elif message.isType("stop film"):
            if self.feed_controller is not None:
                message.addResponse(halMessage.HalMessageResponse(source = self.module_name,
                                                                  data = {"parameters" : self.feed_controller.getParameters()}))

