import sys
import xml.etree.ElementTree as ET  # https://docs.python.org/2/library/xml.etree.elementtree.html
from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QDoubleValidator

class QHLine(QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)

class PKPDParams(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()

        # self.current_param = None
        self.xml_root = None
        label_width = 110
        domain_value_width = 100
        value_width = 60
        label_height = 20
        units_width = 70

        self.scroll = QScrollArea()  # might contain centralWidget

        self.pkpd_params = QWidget()
        self.vbox = QVBoxLayout()
        self.vbox.addStretch(0)
        
        #============ Drug-specific Parameters ================================
        label = QLabel("Drug-specific Parameters")
        label.setFixedHeight(label_height)
        label.setStyleSheet("background-color: orange")
        label.setAlignment(QtCore.Qt.AlignCenter)
        self.vbox.addWidget(label)
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Dosing intervals")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.dose_interval = QLineEdit()
        self.dose_interval.setFixedWidth(label_width)
        self.dose_interval.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.dose_interval)

        units = QLabel("min")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Boundary concentration after dosing")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.dirichlet_node_value_on_dose = QLineEdit()
        self.dirichlet_node_value_on_dose.setFixedWidth(label_width)
        self.dirichlet_node_value_on_dose.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.dirichlet_node_value_on_dose)

        units = QLabel("mmHg")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Decay rate of Dirichlet nodes")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.dirichlet_decay_rate = QLineEdit()
        self.dirichlet_decay_rate.setFixedWidth(label_width)
        self.dirichlet_decay_rate.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.dirichlet_decay_rate)

        units = QLabel("1/min")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Maximum number of doses")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.max_number_doses = QLineEdit()
        self.max_number_doses.setFixedWidth(label_width)
        self.max_number_doses.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.max_number_doses)

        units = QLabel("dimensionless")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)

        #======================================================================
        self.pkpd_params.setLayout(self.vbox)
        self.scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll.setWidgetResizable(True)

        self.scroll.setWidget(self.pkpd_params) # self.config_params = QWidget()
        self.layout = QVBoxLayout(self)
        self.layout.addWidget(self.scroll)
     
    def clear_gui(self):
        pass

    def fill_gui(self):
        print("\n\n------------  pkpd_params_tab: fill_gui --------------")
        uep_pkpd_params = self.xml_root.find(".//pkpd_parameters")
        # custom_data_path = ".//cell_definition[" + str(self.idx_current_cell_def) + "]//custom_data//"
        print('uep_pkpd_params=',uep_pkpd_params)
        if uep_pkpd_params:
            self.dose_interval.setText(self.xml_root.find(".//pkpd_parameters//dose_interval").text)
            self.dirichlet_node_value_on_dose.setText(self.xml_root.find(".//pkpd_parameters//dirichlet_node_value_on_dose").text)
            self.dirichlet_decay_rate.setText(self.xml_root.find(".//pkpd_parameters//dirichlet_decay_rate").text)
            self.max_number_doses.setText(self.xml_root.find(".//pkpd_parameters//max_number_doses").text)
            
    def fill_xml(self):
        print("----------- pkpd_params_tab.py: fill_xml(): ----------")
        self.xml_root.find(".//pkpd_parameters//dose_interval").text = self.dose_interval.text()
        self.xml_root.find(".//pkpd_parameters//dirichlet_node_value_on_dose").text = self.dirichlet_node_value_on_dose.text()
        self.xml_root.find(".//pkpd_parameters//dirichlet_decay_rate").text = self.dirichlet_decay_rate.text()
        self.xml_root.find(".//pkpd_parameters//max_number_doses").text = self.max_number_doses.text()
