# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'qtdesigner/timelapse.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(456, 186)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Dialog.sizePolicy().hasHeightForWidth())
        Dialog.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.channelsBox = QtWidgets.QGroupBox(Dialog)
        self.channelsBox.setObjectName("channelsBox")
        self.gridWidget = QtWidgets.QWidget(self.channelsBox)
        self.gridWidget.setGeometry(QtCore.QRect(10, 30, 511, 41))
        self.gridWidget.setObjectName("gridWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridWidget)
        self.gridLayout.setContentsMargins(5, 12, 5, 20)
        self.gridLayout.setObjectName("gridLayout")
        self.verticalLayout.addWidget(self.channelsBox)
        self.frame = QtWidgets.QFrame(Dialog)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.sliderInterval = QtWidgets.QSlider(self.frame)
        self.sliderInterval.setGeometry(QtCore.QRect(10, 10, 160, 22))
        self.sliderInterval.setMaximum(999)
        self.sliderInterval.setOrientation(QtCore.Qt.Horizontal)
        self.sliderInterval.setObjectName("sliderInterval")
        self.comboIntervalType = QtWidgets.QComboBox(self.frame)
        self.comboIntervalType.setGeometry(QtCore.QRect(280, 6, 91, 30))
        self.comboIntervalType.setObjectName("comboIntervalType")
        self.comboIntervalType.addItem("")
        self.comboIntervalType.addItem("")
        self.comboIntervalType.addItem("")
        self.comboIntervalType.addItem("")
        self.textInterval = QtWidgets.QLineEdit(self.frame)
        self.textInterval.setGeometry(QtCore.QRect(180, 10, 61, 21))
        self.textInterval.setObjectName("textInterval")
        self.verticalLayout.addWidget(self.frame)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem)
        self.okButton = QtWidgets.QPushButton(Dialog)
        self.okButton.setObjectName("okButton")
        self.horizontalLayout_3.addWidget(self.okButton)
        self.verticalLayout.addLayout(self.horizontalLayout_3)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Timelapse setup"))
        self.channelsBox.setTitle(_translate("Dialog", "Channels"))
        self.comboIntervalType.setItemText(0, _translate("Dialog", "ms"))
        self.comboIntervalType.setItemText(1, _translate("Dialog", "s"))
        self.comboIntervalType.setItemText(2, _translate("Dialog", "min"))
        self.comboIntervalType.setItemText(3, _translate("Dialog", "h"))
        self.textInterval.setText(_translate("Dialog", "0"))
        self.okButton.setText(_translate("Dialog", "Ok"))

