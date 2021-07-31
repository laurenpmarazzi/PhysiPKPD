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
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Drug concentration in systemic circulation")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.systemic_circulation_increase_on_dose = QLineEdit()
        self.systemic_circulation_increase_on_dose.setFixedWidth(label_width)
        self.systemic_circulation_increase_on_dose.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.systemic_circulation_increase_on_dose)

        units = QLabel("mmHg")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Drug decay rate in systemic circulation")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.systemic_circulation_elimination_rate = QLineEdit()
        self.systemic_circulation_elimination_rate.setFixedWidth(label_width)
        self.systemic_circulation_elimination_rate.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.systemic_circulation_elimination_rate)

        units = QLabel("1/min")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Volume ratio of systemic circulation to periphery")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.systemic_circulation_to_periphery_volume_ratio = QLineEdit()
        self.systemic_circulation_to_periphery_volume_ratio.setFixedWidth(label_width)
        self.systemic_circulation_to_periphery_volume_ratio.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.systemic_circulation_to_periphery_volume_ratio)

        units = QLabel("dimensionless")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Biotransport ratio")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.biot_number = QLineEdit()
        self.biot_number.setFixedWidth(label_width)
        self.biot_number.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.biot_number)

        units = QLabel("dimensionless")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Drug flux across capillaries")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.drug_flux_across_capillaries = QLineEdit()
        self.drug_flux_across_capillaries.setFixedWidth(label_width)
        self.drug_flux_across_capillaries.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.drug_flux_across_capillaries)

        units = QLabel("1/min")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)
        
        #----------
        hbox = QHBoxLayout()
        self.set_first_dose_time = QCheckBox("User sets first dose time")
        # self.motility_enabled.setAlignment(QtCore.Qt.AlignRight)
        # label.setFixedWidth(label_width)
        hbox.addWidget(self.set_first_dose_time)
        self.vbox.addLayout(hbox)
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Time of first dose")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.first_dose_time = QLineEdit()
        self.first_dose_time.setFixedWidth(label_width)
        self.first_dose_time.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.first_dose_time)

        units = QLabel("min")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Confluence (between 0 to 1) to start therapy")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.confluence_condition = QLineEdit()
        self.confluence_condition.setFixedWidth(label_width)
        self.confluence_condition.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.confluence_condition)

        units = QLabel("dimensionless")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)

        #============ Apoptosis Parameters ================================
        label = QLabel("Apoptosis Parameters")
        label.setFixedHeight(label_height)
        label.setStyleSheet("background-color: orange")
        label.setAlignment(QtCore.Qt.AlignCenter)
        self.vbox.addWidget(label)
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Apoptosis EC50")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.EC_50 = QLineEdit()
        self.EC_50.setFixedWidth(label_width)
        self.EC_50.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.EC_50)

        units = QLabel("dimensionless")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)

        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Apoptosis hill coefficient")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.Hill_power = QLineEdit()
        self.Hill_power.setFixedWidth(label_width)
        self.Hill_power.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.Hill_power)

        units = QLabel("dimensionless")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Maximum increase to apoptosis rate with drug")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.max_increase_to_apoptosis = QLineEdit()
        self.max_increase_to_apoptosis.setFixedWidth(label_width)
        self.max_increase_to_apoptosis.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.max_increase_to_apoptosis)

        units = QLabel("1/min")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)
        
        #--------------------
        hbox = QHBoxLayout()
        label = QLabel("Use Hill function for apoptosis damage (set 1 for true)")
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.use_AUC_into_hill = QLineEdit()
        self.use_AUC_into_hill.setFixedWidth(label_width)
        self.use_AUC_into_hill.setValidator(QtGui.QDoubleValidator())
        hbox.addWidget(self.use_AUC_into_hill)

        units = QLabel("dimensionless")
        units.setAlignment(QtCore.Qt.AlignLeft)
        units.setFixedWidth(units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)
        
        #============ Parameters for Mechanism of Action ================================
        label = QLabel("Parameters for Mechanism of Action")
        label.setFixedHeight(label_height)
        label.setStyleSheet("background-color: orange")
        label.setAlignment(QtCore.Qt.AlignCenter)
        self.vbox.addWidget(label)
        
        #----------
        hbox = QHBoxLayout()
        self.moa_apoptosis = QCheckBox("MOA - Apoptotic signaling")
        # self.motility_enabled.setAlignment(QtCore.Qt.AlignRight)
        # label.setFixedWidth(label_width)
        hbox.addWidget(self.moa_apoptosis)
        self.vbox.addLayout(hbox)
        
        #----------
        hbox = QHBoxLayout()
        self.moa_proliferation = QCheckBox("MOA - Proliferation arrest")
        # self.motility_enabled.setAlignment(QtCore.Qt.AlignRight)
        # label.setFixedWidth(label_width)
        hbox.addWidget(self.moa_proliferation)
        self.vbox.addLayout(hbox)
        
        #----------
        hbox = QHBoxLayout()
        self.moa_necrosis = QCheckBox("MOA - Necrotic signaling")
        # self.motility_enabled.setAlignment(QtCore.Qt.AlignRight)
        # label.setFixedWidth(label_width)
        hbox.addWidget(self.moa_necrosis)
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
            self.systemic_circulation_increase_on_dose.setText(self.xml_root.find(".//pkpd_parameters//systemic_circulation_increase_on_dose").text)
            self.systemic_circulation_elimination_rate.setText(self.xml_root.find(".//pkpd_parameters//systemic_circulation_elimination_rate").text)
            self.systemic_circulation_to_periphery_volume_ratio.setText(self.xml_root.find(".//pkpd_parameters//systemic_circulation_to_periphery_volume_ratio").text)
            self.biot_number.setText(self.xml_root.find(".//pkpd_parameters//biot_number").text)
            self.drug_flux_across_capillaries.setText(self.xml_root.find(".//pkpd_parameters//drug_flux_across_capillaries").text)
            self.EC_50.setText(self.xml_root.find(".//pkpd_parameters//EC_50").text)
            self.Hill_power.setText(self.xml_root.find(".//pkpd_parameters//Hill_power").text)
            self.max_increase_to_apoptosis.setText(self.xml_root.find(".//pkpd_parameters//max_increase_to_apoptosis").text)
            self.use_AUC_into_hill.setText(self.xml_root.find(".//pkpd_parameters//use_AUC_into_hill").text)
            self.first_dose_time.setText(self.xml_root.find(".//pkpd_parameters//first_dose_time").text)
            self.confluence_condition.setText(self.xml_root.find(".//pkpd_parameters//confluence_condition").text)
            
            # Set boolean
            if self.xml_root.find(".//pkpd_parameters//moa_apoptosis").text.lower() == 'true':
                self.moa_apoptosis.setChecked(True)
            else:
                self.moa_apoptosis.setChecked(False)
                
            if self.xml_root.find(".//pkpd_parameters//moa_proliferation").text.lower() == 'true':
                self.moa_proliferation.setChecked(True)
            else:
                self.moa_proliferation.setChecked(False)
                
            if self.xml_root.find(".//pkpd_parameters//moa_necrosis").text.lower() == 'true':
                self.moa_necrosis.setChecked(True)
            else:
                self.moa_necrosis.setChecked(False)
                
            if self.xml_root.find(".//pkpd_parameters//set_first_dose_time").text.lower() == 'true':
                self.set_first_dose_time.setChecked(True)
            else:
                self.set_first_dose_time.setChecked(False)

    def fill_xml(self):
        print("----------- pkpd_params_tab.py: fill_xml(): ----------")
        self.xml_root.find(".//pkpd_parameters//dose_interval").text = self.dose_interval.text()
        self.xml_root.find(".//pkpd_parameters//dirichlet_node_value_on_dose").text = self.dirichlet_node_value_on_dose.text()
        self.xml_root.find(".//pkpd_parameters//dirichlet_decay_rate").text = self.dirichlet_decay_rate.text()
        self.xml_root.find(".//pkpd_parameters//max_number_doses").text = self.max_number_doses.text()
        
        self.xml_root.find(".//pkpd_parameters//systemic_circulation_increase_on_dose").text = self.systemic_circulation_increase_on_dose.text()
        self.xml_root.find(".//pkpd_parameters//systemic_circulation_elimination_rate").text = self.systemic_circulation_elimination_rate.text()
        self.xml_root.find(".//pkpd_parameters//systemic_circulation_to_periphery_volume_ratio").text = self.systemic_circulation_to_periphery_volume_ratio.text()
        self.xml_root.find(".//pkpd_parameters//biot_number").text = self.biot_number.text()
        self.xml_root.find(".//pkpd_parameters//drug_flux_across_capillaries").text = self.drug_flux_across_capillaries.text()
        self.xml_root.find(".//pkpd_parameters//EC_50").text = self.EC_50.text()
        self.xml_root.find(".//pkpd_parameters//Hill_power").text = self.Hill_power.text()
        self.xml_root.find(".//pkpd_parameters//max_increase_to_apoptosis").text = self.max_increase_to_apoptosis.text()
        
        self.xml_root.find(".//pkpd_parameters//use_AUC_into_hill").text = self.use_AUC_into_hill.text()
        self.xml_root.find(".//pkpd_parameters//first_dose_time").text = self.first_dose_time.text()
        self.xml_root.find(".//pkpd_parameters//confluence_condition").text = self.confluence_condition.text()
        
        bval = "false"
        if self.moa_apoptosis.isChecked():
            bval = "true"
        self.xml_root.find(".//pkpd_parameters//moa_apoptosis").text = bval
        
        bval = "false"
        if self.moa_proliferation.isChecked():
            bval = "true"
        self.xml_root.find(".//pkpd_parameters//moa_proliferation").text = bval
        
        bval = "false"
        if self.moa_necrosis.isChecked():
            bval = "true"
        self.xml_root.find(".//pkpd_parameters//moa_necrosis").text = bval
        
        bval = "false"
        if self.set_first_dose_time.isChecked():
            bval = "true"
        self.xml_root.find(".//pkpd_parameters//set_first_dose_time").text = bval
