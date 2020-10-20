from PyQt5 import QtWidgets, QtGui, uic
from PyQt5.QtCore import QRunnable, QThreadPool, pyqtSlot, QTimer
from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
import numpy as np
import joblib

import sys
import os
import subprocess
import shutil
import xml.etree.ElementTree as ET
from xml.dom import minidom

import classes
import schemes
import fitting
import detection

class Ui(QtWidgets.QMainWindow):

    def __init__(self):
        super(Ui, self).__init__()  # Call the inherited classes __init__ method
        exe = os.path.dirname(__file__)
        uic.loadUi(exe + '/gui.ui', self)  # Load the .ui file
        
        # Buttons actions
        self.toolButton_input_file.clicked.connect(self.browse_input_file_pressed)
        self.toolButton_mean_file.clicked.connect(self.browse_mean_file_pressed)
        self.toolButton_output_folder.clicked.connect(self.browse_output_folder_pressed)
        self.pushButton_run.clicked.connect(self.run_pressed)
        self.actionOpen.triggered.connect(self.open_config_clicked)
        self.actionSave.triggered.connect(self.save_config_clicked)
        self.actionQuit.triggered.connect(self.quit_clicked)

        self.cwd = os.getcwd()
        self.lineEdit_output_folder.setText(self.cwd)

        # Detection Frame
        self.fig1 = Figure()
        self.canvas1 = FigureCanvas(self.fig1)
        self.toolbar1 = NavigationToolbar(self.canvas1, self)
        self.gridLayout_detection.addWidget(self.toolbar1)
        self.gridLayout_detection.addWidget(self.canvas1)

        # Fitting Frame
        self.fig2 = Figure()
        self.canvas2 = FigureCanvas(self.fig2)
        self.toolbar2 = NavigationToolbar(self.canvas2, self)
        self.gridLayout_fitting.addWidget(self.toolbar2)
        self.gridLayout_fitting.addWidget(self.canvas2)
        
        self.show()

    def quit_clicked(self):
        sys.exit(0)

    def open_config_clicked(self):
        filename = self.open_config("xml")
        self.open_config_xml(filename)

    def open_config_xml(self, xml_file):
        tree = ET.parse(xml_file)
        for attr in tree.iter():
            if attr.tag == 'VortexFitting':
                continue
            while True:
                obj = self.findChild(QtWidgets.QLineEdit, attr.tag)
                if obj:
                    obj.setText(attr.text)
                    break
                obj = self.findChild(QtWidgets.QComboBox, attr.tag)
                if obj:
                    obj.setCurrentText(attr.text)
                    break
                obj = self.findChild(QtWidgets.QCheckBox, attr.tag)
                if obj:
                    obj.setChecked(True)
                    break
                obj = self.findChild(QtWidgets.QGroupBox, attr.tag)
                if obj:
                    if attr.text == "True":
                        obj.setChecked(True)
                    break
                else:
                    break

    def save_config_clicked(self):

        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save File")
        data = ET.Element('VortexFitting')

        for attr, value in self.__dict__.items():
            if "lineEdit_case_name" in attr:
                continue
            if value.__class__.__name__ == 'QLineEdit':
                ET.SubElement(data, str(attr)).text = value.text()
            elif value.__class__.__name__ == 'QComboBox':
                ET.SubElement(data, str(attr)).text = value.currentText()
            elif value.__class__.__name__ == 'QCheckBox':
                ET.SubElement(data, str(attr)).text = str(value.isChecked())
            elif value.__class__.__name__ == "QGroupBox":
                ET.SubElement(data, str(attr)).text = str(value.isChecked())

        xmlstr = minidom.parseString(ET.tostring(data)).toprettyxml(indent="   ")
        with open(filename, "w") as f:
            f.write(xmlstr)
    
    
    def missing_file_popup(self, filename):
        print("File {} not found!".format(filename))
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Critical)
        msg.setWindowTitle("Error")
        msg.setText("File {} not found!".format(filename))
        msg.exec_()

    def missing_option(self, option):
        print("Missing {} option!".format(option))
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Critical)
        msg.setWindowTitle("Error")
        msg.setText("Missing {} option!".format(option))
        msg.exec_()

    def open_config(self, ext):
        options = QtWidgets.QFileDialog.Options()
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open", "", "{}(*.{})".format(ext,ext), options=options)
        return filename

    def browse_input_file_pressed(self):
        options = QtWidgets.QFileDialog.Options()
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open", options=options)
        self.lineEdit_input_file.setText(filename)

    def browse_mean_file_pressed(self):
        options = QtWidgets.QFileDialog.Options()
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open", options=options)
        self.lineEdit_mean_file.setText(filename)

    def browse_output_folder_pressed(self):
        foldername = QtWidgets.QFileDialog.getExistingDirectory(self, "Select Output Folder")
        if foldername:
            self.lineEdit_output_folder.setText(foldername)
            os.chdir(foldername)        


    def add_fig_detection(self, vortices_counterclockwise, vortices_clockwise, detection_field, flip_axis):
        ax = self.fig1.add_subplot(111)
        ax.clear()

        if flip_axis:
            detection_field = detection_field.T  # transpose the detection field
            ax.scatter(vortices_counterclockwise[0], vortices_counterclockwise[1],
                        edgecolor='green', facecolor='green', label='counterclockwise')
            ax.scatter(vortices_clockwise[0], vortices_clockwise[1],
                        edgecolor='yellow', facecolor='yellow', label='clockwise')
        else:
            ax.scatter(vortices_counterclockwise[1], vortices_counterclockwise[0],
                        edgecolor='green', facecolor='green', label='counterclockwise')
            ax.scatter(vortices_clockwise[1], vortices_clockwise[0],
                        edgecolor='yellow', facecolor='yellow', label='clockwise')

        ax.set_title('Detected possible vortices')
        # plt.contourf(field, cmap="Greys_r")

        ax.imshow(detection_field, origin='lower', cmap="Greys_r")
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        # plt.legend()
        # plt.imshow(field, cmap="Greys_r",origin="lower")

        #plt.show()
        self.canvas1.draw()

    def add_fig_fitting(self):
        ax = self.fig2.add_subplot(111)
        ax.clear()
        x = [0.1, 0.2, 0.3]
        y = [5, 4, 3]
        ax.plot(x, y)
        ax.set_xlabel('x')
        ax.set_ylabel('y')

        self.canvas2.draw()


    def run_pressed(self):

        vfield = classes.VelocityField(self.lineEdit_input_file.text(),
                                       self.lineEdit_timestep.text(),
                                       self.lineEdit_mean_file.text(),
                                       self.comboBox_filetype.currentText())

        # ---- DIFFERENCE APPROXIMATION ----#
        if self.comboBox_scheme.currentText() == 'Forth-Order':
            vfield.derivative = schemes.fourth_order_diff(vfield)
        elif self.comboBox_scheme.currentText() == 'Second-Order':
            vfield.derivative = schemes.second_order_diff(vfield)
        elif self.comboBox_scheme.currentText() == 'Least-square filter':
            vfield.derivative = schemes.least_square_diff(vfield)
        else:
            print('No scheme', args.scheme, 'found. Exiting!')
            missing_option(comboBox_scheme.currentText())

        # ---- VORTICITY ----#
        vorticity = vfield.derivative['dvdx'] - vfield.derivative['dudy']

        # ---- METHOD FOR DETECTION OF VORTICES ----#
        detection_field = []
        if self.comboBox_criterion.currentText() == 'Q':
            detection_field = detection.calc_q_criterion(vfield)
        elif self.comboBox_criterion.currentText() == 'swirling':
            detection_field = detection.calc_swirling(vfield)
        elif self.comboBox_criterion.currentText() == 'delta':
            detection_field = detection.calc_delta_criterion(vfield)

        if self.checkBox_normalization.isChecked():
            detection_field = fitting.normalize(detection_field, vfield.normalization_direction)  # normalization

        # ---- PEAK DETECTION ----#
        print('Threshold=', self.lineEdit_threshold.text(),
              ', box size=', self.lineEdit_boxsize.text())

        peaks = fitting.find_peaks(detection_field, float(self.lineEdit_threshold.text()),
                                   int(self.lineEdit_boxsize.text()))

        print('Vortices found: ', len(peaks[0]))
        # ---- PEAKS DIRECTION OF ROTATION ----#
        vortices_counterclockwise, vortices_clockwise = fitting.direction_rotation(vorticity, peaks)

        # ---- MODEL FITTING ----#
        vortices = list()
        if self.groupBox_fitting.isChecked():
            vortices = fitting.get_vortices(vfield, peaks, vorticity,
                                            float(self.lineEdit_max_radius.text()),
                                            float(self.lineEdit_threshold.text()))
            print('---- Accepted vortices ----')
            print(len(vortices))
        else:
            print('No fitting')

        self.add_fig_detection(vortices_counterclockwise, vortices_clockwise, detection_field, self.checkBox_flip_axis.isChecked())
        self.add_fig_fitting()
        # Change to calibration preview
        #case.tabWidget_right.setCurrentIndex(0)


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)  # Create an instance of QtWidgets.QApplication
    window = Ui()  # Create an instance of our class

    timer = QTimer()
    timer.timeout.connect(lambda: None)
    timer.start(100)

    app.exec_()  # Start the application
