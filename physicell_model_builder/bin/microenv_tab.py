"""
Authors:
Randy Heiland (heiland@iu.edu)
Adam Morrow, Grant Waldrow, Drew Willis, Kim Crevecoeur
Dr. Paul Macklin (macklinp@iu.edu)

--- Versions ---
0.1 - initial version
"""

import sys
import copy
import xml.etree.ElementTree as ET  # https://docs.python.org/2/library/xml.etree.elementtree.html
# from ElementTree_pretty import prettify

from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QDoubleValidator #, QTreeWidgetItemIterator

class QHLine(QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)

class SubstrateDef(QWidget):
    def __init__(self):
        super().__init__()
        # global self.microenv_params

        self.param_d = {}  # a dict of dicts - rwh/todo, used anymore?
        # self.substrate = {}
        self.current_substrate = None
        self.xml_root = None
        self.celldef_tab = None
        self.new_substrate_count = 1

        # self.stacked_w = QStackedWidget()
        # self.stack_w = []
        # self.stack_w.append(QStackedWidget())
        # self.stacked_w.addWidget(self.stack_w[0])

        #---------------
        # self.cell_defs = CellDefInstances()
        self.microenv_hbox = QHBoxLayout()

        self.splitter = QSplitter()

        tree_widget_width = 160
        tree_widget_height = 400

        self.tree = QTreeWidget() # tree is overkill; list would suffice; Meh.
        # self.tree.itemDoubleClicked.connect(self.treeitem_edit_cb)
        # self.tree.setStyleSheet("background-color: lightgray")
        self.tree.setFixedWidth(tree_widget_width)
        self.tree.setFixedHeight(tree_widget_height)
        # self.tree.currentChanged(self.tree_item_changed_cb)
        self.tree.itemClicked.connect(self.tree_item_clicked_cb)
        # self.tree.itemSelectionChanged.connect(self.tree_item_changed_cb2)
        # self.tree.itemDoubleClicked.connect(self.tree_item_changed_cb2)
        # self.tree.itemSelectionChanged.connect(self.tree_item_changed_cb2)
        self.tree.itemChanged.connect(self.tree_item_changed_cb)   # rename a substrate
        # self.tree.selectionChanged.connect(self.tree_item_sel_changed_cb)
        # self.tree.currentChanged.connect(self.tree_item_sel_changed_cb)
        # self.tree.itemSelectionChanged()
        # self.tree.setColumnCount(1)

        # self.tree.setCurrentItem(0)  # rwh/TODO

        header = QTreeWidgetItem(["---  Substrate ---"])
        self.tree.setHeaderItem(header)

        # cellname = QTreeWidgetItem(["virus"])
        # self.tree.insertTopLevelItem(0,cellname)

        # cellname = QTreeWidgetItem(["interferon"])
        # self.tree.insertTopLevelItem(1,cellname)

        #-------------------------
        # self.name_list = QListWidget() # tree is overkill; list would suffice; meh.

        self.microenv_hbox.addWidget(self.tree)

        # self.microenv_hbox.addWidget(self.name_list)


        self.scroll_cell_def_tree = QScrollArea()
        self.scroll_cell_def_tree.setWidget(self.tree)
        # self.scroll_cell_def_tree.setWidget(self.name_list)

        # splitter.addWidget(self.tree)
        self.splitter.addWidget(self.scroll_cell_def_tree)

        #-------------------------------------------
        # self.tab = QWidget()
        # self.tabs.resize(200,5)
        
        #-------------------------------------------
        self.label_width = 150
        self.units_width = 80
        self.tab_widget =  QTabWidget()
        # self.scroll = QScrollArea()
        self.scroll_area = QScrollArea()
        self.splitter.addWidget(self.scroll_area)
        # self.microenv_hbox.addWidget(self.scroll_area)

        self.microenv_params = QWidget()
        self.vbox = QVBoxLayout()
        self.vbox.addStretch(0)

        # self.microenv_hbox.addWidget(self.)

        #------------------
        controls_hbox = QHBoxLayout()
        self.new_button = QPushButton("New")
        self.new_button.clicked.connect(self.new_substrate)
        controls_hbox.addWidget(self.new_button)

        self.copy_button = QPushButton("Copy")
        self.copy_button.clicked.connect(self.copy_substrate)
        controls_hbox.addWidget(self.copy_button)

        self.delete_button = QPushButton("Delete")
        self.delete_button.clicked.connect(self.delete_substrate)
        controls_hbox.addWidget(self.delete_button)

        self.main_tab = QWidget()
        self.secretion_tab = QWidget()
        self.scroll_params = QScrollArea()
        #self.splitter.addWidget(self.scroll_params)
        self.tab_widget.addTab(self.create_main_tab(),"main")
        self.tab_widget.addTab(self.create_secretion_tab(),"Secretion")
        self.microenv_tabs_layout = QGridLayout()
        self.microenv_tabs_layout.addWidget(self.tab_widget, 0,0,1,1) # w, row, column, rowspan, colspan
        # #------------------
        # hbox = QHBoxLayout()
        # label = QLabel("diffusion coefficient")
        # label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # hbox.addWidget(label)

        # self.diffusion_coef = QLineEdit()
        # self.diffusion_coef.setValidator(QtGui.QDoubleValidator())
        # self.diffusion_coef.textChanged.connect(self.diffusion_coef_changed)
        # # self.diffusion_coef.enter.connect(self.save_xml)
        # hbox.addWidget(self.diffusion_coef)

        # units = QLabel("micron^2/min")
        # units.setFixedWidth(self.units_width)
        # hbox.addWidget(units)
        # self.vbox.addLayout(hbox)

        # #----------
        # hbox = QHBoxLayout()
        # label = QLabel("decay rate")
        # label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # hbox.addWidget(label)

        # self.decay_rate = QLineEdit()
        # self.decay_rate.setValidator(QtGui.QDoubleValidator())
        # self.decay_rate.textChanged.connect(self.decay_rate_changed)
        # # self.decay_rate.enter.connect(self.save_xml)
        # hbox.addWidget(self.decay_rate)

        # units = QLabel("1/min")
        # units.setFixedWidth(self.units_width)
        # hbox.addWidget(units)
        # self.vbox.addLayout(hbox)

        # #----------
        # hbox = QHBoxLayout()
        # label = QLabel("initial condition")
        # label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # hbox.addWidget(label)

        # self.init_cond = QLineEdit()
        # self.init_cond.setValidator(QtGui.QDoubleValidator())
        # self.init_cond.textChanged.connect(self.init_cond_changed)
        # # self.init_cond.enter.connect(self.save_xml)
        # hbox.addWidget(self.init_cond)

        # units = QLabel("mmol")
        # units.setFixedWidth(self.units_width)
        # hbox.addWidget(units)
        # self.vbox.addLayout(hbox)
        # #----------

        # hbox = QHBoxLayout()
        # label = QLabel("Dirichlet BC")
        # label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # hbox.addWidget(label)

        # self.dirichlet_bc = QLineEdit()
        # self.dirichlet_bc.setValidator(QtGui.QDoubleValidator())
        # self.dirichlet_bc.textChanged.connect(self.dirichlet_bc_changed)
        # # self.bdy_cond.enter.connect(self.save_xml)
        # hbox.addWidget(self.dirichlet_bc)

        # units = QLabel("mmol")
        # units.setFixedWidth(self.units_width)
        # hbox.addWidget(units)

        # self.dirichlet_bc_enabled = QCheckBox("on")
        # self.dirichlet_bc_enabled.stateChanged.connect(self.dirichlet_toggle_cb)
        # # self.motility_enabled.setAlignment(QtCore.Qt.AlignRight)
        # # label.setFixedWidth(self.label_width)
        # hbox.addWidget(self.dirichlet_bc_enabled)

        # self.vbox.addLayout(hbox)

        # #--------------------------

        # dirichlet_options_bdy = QLabel("Dirichlet options per boundary:")
        # # units.setFixedWidth(self.units_width)
        # self.vbox.addWidget(dirichlet_options_bdy)

        # #----
        # hbox = QHBoxLayout()
        # label = QLabel("xmin:")
        # label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # hbox.addWidget(label)

        # self.dirichlet_xmin = QLineEdit()
        # self.dirichlet_xmin.setValidator(QtGui.QDoubleValidator())
        # self.dirichlet_xmin.textChanged.connect(self.dirichlet_xmin_changed)
        # hbox.addWidget(self.dirichlet_xmin)

        # self.enable_xmin = QCheckBox("on")
        # self.enable_xmin.stateChanged.connect(self.enable_xmin_cb)
        # # self.motility_enabled.setAlignment(QtCore.Qt.AlignRight)
        # # label.setFixedWidth(self.label_width)
        # hbox.addWidget(self.enable_xmin)
        # self.vbox.addLayout(hbox)
        # #----
        # hbox = QHBoxLayout()
        # label = QLabel("xmax:")
        # label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # hbox.addWidget(label)

        # self.dirichlet_xmax = QLineEdit()
        # self.dirichlet_xmax.setValidator(QtGui.QDoubleValidator())
        # self.dirichlet_xmax.textChanged.connect(self.dirichlet_xmax_changed)
        # hbox.addWidget(self.dirichlet_xmax)

        # self.enable_xmax = QCheckBox("on")
        # self.enable_xmax.stateChanged.connect(self.enable_xmax_cb)
        # hbox.addWidget(self.enable_xmax)
        # self.vbox.addLayout(hbox)
        # #---------
        # hbox = QHBoxLayout()
        # label = QLabel("ymin:")
        # label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # hbox.addWidget(label)

        # self.dirichlet_ymin = QLineEdit()
        # self.dirichlet_ymin.setValidator(QtGui.QDoubleValidator())
        # self.dirichlet_ymin.textChanged.connect(self.dirichlet_ymin_changed)
        # hbox.addWidget(self.dirichlet_ymin)

        # self.enable_ymin = QCheckBox("on")
        # self.enable_ymin.stateChanged.connect(self.enable_ymin_cb)
        # # self.motility_enabled.setAlignment(QtCore.Qt.AlignRight)
        # # label.setFixedWidth(self.label_width)
        # hbox.addWidget(self.enable_ymin)
        # self.vbox.addLayout(hbox)
        # #----
        # hbox = QHBoxLayout()
        # label = QLabel("ymax:")
        # label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # hbox.addWidget(label)

        # self.dirichlet_ymax = QLineEdit()
        # self.dirichlet_ymax.setValidator(QtGui.QDoubleValidator())
        # self.dirichlet_ymax.textChanged.connect(self.dirichlet_ymax_changed)
        # hbox.addWidget(self.dirichlet_ymax)

        # self.enable_ymax = QCheckBox("on")
        # self.enable_ymax.stateChanged.connect(self.enable_ymax_cb)
        # hbox.addWidget(self.enable_ymax)
        # self.vbox.addLayout(hbox)
        # #---------
        # hbox = QHBoxLayout()
        # label = QLabel("zmin:")
        # label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # hbox.addWidget(label)

        # self.dirichlet_zmin = QLineEdit()
        # self.dirichlet_zmin.setValidator(QtGui.QDoubleValidator())
        # self.dirichlet_zmin.textChanged.connect(self.dirichlet_zmin_changed)
        # hbox.addWidget(self.dirichlet_zmin)

        # self.enable_zmin = QCheckBox("on")
        # self.enable_zmin.stateChanged.connect(self.enable_zmin_cb)
        # # self.motility_enabled.setAlignment(QtCore.Qt.AlignRight)
        # # label.setFixedWidth(self.label_width)
        # hbox.addWidget(self.enable_zmin)
        # self.vbox.addLayout(hbox)
        # #----
        # hbox = QHBoxLayout()
        # label = QLabel("zmax:")
        # label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # hbox.addWidget(label)

        # self.dirichlet_zmax = QLineEdit()
        # self.dirichlet_zmax.setValidator(QtGui.QDoubleValidator())
        # self.dirichlet_zmax.textChanged.connect(self.dirichlet_zmax_changed)
        # hbox.addWidget(self.dirichlet_zmax)

        # self.enable_zmax = QCheckBox("on")
        # self.enable_zmax.stateChanged.connect(self.enable_zmax_cb)
        # hbox.addWidget(self.enable_zmax)
        # self.vbox.addLayout(hbox)

        # #-------------
        # # Toggles for overall microenv (all substrates)
        # self.vbox.addWidget(QHLine())

        # hbox = QHBoxLayout()
        # hbox.addWidget(QLabel("For all substrates: "))

        # self.gradients = QCheckBox("calculate gradients")
        # self.gradients.stateChanged.connect(self.gradients_cb)
        # hbox.addWidget(self.gradients)

        # self.track_in_agents = QCheckBox("track in agents")
        # self.track_in_agents.stateChanged.connect(self.track_in_agents_cb)
        # hbox.addWidget(self.track_in_agents)
        # self.vbox.addLayout(hbox)

        # self.vbox.addStretch()

        # self.microenv_params.setLayout(self.vbox)

        # self.scroll_area.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        # self.scroll_area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        # self.scroll_area.setWidgetResizable(True)
        # self.scroll_area.setWidget(self.microenv_params)

        # self.layout = QVBoxLayout(self)

        self.layout.addLayout(controls_hbox)

        self.layout.addWidget(self.splitter)


    def diffusion_coef_changed(self, text):
        print("Text: %s", text)
        self.param_d[self.current_substrate]["diffusion_coef"] = text
        # log.info("diffusion_coef changed: %s", text)

    def decay_rate_changed(self, text):
        self.param_d[self.current_substrate]["decay_rate"] = text
    def init_cond_changed(self, text):
        self.param_d[self.current_substrate]["init_cond"] = text
    def dirichlet_bc_changed(self, text):
        self.param_d[self.current_substrate]["dirichlet_bc"] = text

    def dirichlet_toggle_cb(self):
        print("dirichlet_toggle_cb()")
        self.param_d[self.current_substrate]["dirichlet_enabled"] = self.dirichlet_bc_enabled.isChecked()

    # global to all substrates
    def gradients_cb(self):
        # self.param_d[self.current_substrate]["gradients"] = self.gradients.isChecked()
        self.param_d["gradients"] = self.gradients.isChecked()
    def track_in_agents_cb(self):
        self.param_d["track_in_agents"] = self.track_in_agents.isChecked()

    def dirichlet_xmin_changed(self, text):
        # print("\n\n------> def dirichlet_xmin_changed(self, text)  called!!!")
        self.param_d[self.current_substrate]["dirichlet_xmin"] = text
    def dirichlet_xmax_changed(self, text):
        self.param_d[self.current_substrate]["dirichlet_xmax"] = text
    def dirichlet_ymin_changed(self, text):
        self.param_d[self.current_substrate]["dirichlet_ymin"] = text
    def dirichlet_ymax_changed(self, text):
        self.param_d[self.current_substrate]["dirichlet_ymax"] = text
    def dirichlet_zmin_changed(self, text):
        self.param_d[self.current_substrate]["dirichlet_zmin"] = text
    def dirichlet_zmax_changed(self, text):
        self.param_d[self.current_substrate]["dirichlet_zmax"] = text

    def enable_xmin_cb(self):
        self.param_d[self.current_substrate]["enable_xmin"] = self.enable_xmin.isChecked()
    def enable_xmax_cb(self):
        self.param_d[self.current_substrate]["enable_xmax"] = self.enable_xmax.isChecked()
    def enable_ymin_cb(self):
        self.param_d[self.current_substrate]["enable_ymin"] = self.enable_ymin.isChecked()
    def enable_ymax_cb(self):
        self.param_d[self.current_substrate]["enable_ymax"] = self.enable_ymax.isChecked()
    def enable_zmin_cb(self):
        self.param_d[self.current_substrate]["enable_zmin"] = self.enable_zmin.isChecked()
    def enable_zmax_cb(self):
        self.param_d[self.current_substrate]["enable_zmax"] = self.enable_zmax.isChecked()

    #----------------------------------------------------------------------
    # @QtCore.Slot()
    def new_substrate(self):
        # print('------ new_substrate')
        subname = "substrate%02d" % self.new_substrate_count
        # Make a new substrate (that's a copy of the currently selected one)
        # self.param_d[subname] = self.param_d[self.current_substrate].copy()  #rwh - "copy()" is critical

        self.param_d[subname] = copy.deepcopy(self.param_d[self.current_substrate])

        # self.param_d[subname]["name"] = subname
        # for k in self.param_d.keys():
        #     print(" (pre-new vals)===>>> ",k, " : ", self.param_d[k])
        # print()

        # Then "zero out" all entries(?)
        text = "0.0"
        self.param_d[subname]["diffusion_coef"] = text
        self.param_d[subname]["decay_rate"] = text
        self.param_d[subname]["init_cond"] = text
        self.param_d[subname]["dirichlet_bc"] = text
        bval = False
        self.param_d[subname]["dirichlet_enabled"] = bval

        text = ""
        self.param_d[subname]["dirichlet_xmin"] = text
        self.param_d[subname]["dirichlet_xmax"] = text
        self.param_d[subname]["dirichlet_ymin"] = text
        self.param_d[subname]["dirichlet_ymax"] = text
        self.param_d[subname]["dirichlet_zmin"] = text
        self.param_d[subname]["dirichlet_zmax"] = text

        bval = False
        self.param_d[subname]["enable_xmin"] = bval
        self.param_d[subname]["enable_xmax"] = bval
        self.param_d[subname]["enable_ymin"] = bval
        self.param_d[subname]["enable_ymax"] = bval
        self.param_d[subname]["enable_zmin"] = bval
        self.param_d[subname]["enable_zmax"] = bval

        self.param_d["gradients"] = bval
        self.param_d["track_in_agents"] = bval

        # print("\n ----- new dict:")
        # for k in self.param_d.keys():
        #     print(" ===>>> ",k, " : ", self.param_d[k])

        self.new_substrate_count += 1

        self.celldef_tab.add_new_substrate(subname)
        # self.celldef_tab.add_new_substrate_comboboxes(subname)
        # self.param_d[cell_def_name]["secretion"][substrate_name] = {}

        # sval = "0.0"
        # print("cdnames (keys) = ",self.celldef_tab.param_d.keys())
        # for cdname in self.celldef_tab.param_d.keys():  # for all cell defs, initialize secretion params
        #     # self.param_d[cdname]["secretion"][self.current_secretion_substrate]["secretion_rate"] = sval
        #     print('cdname = ',cdname)
        #     print(self.celldef_tab.param_d[cdname]["secretion"])
        #     self.celldef_tab.param_d[cdname]["secretion"][subname]["secretion_rate"] = sval
        #     self.celldef_tab.param_d[cdname]["secretion"][subname]["secretion_target"] = sval
        #     self.celldef_tab.param_d[cdname]["secretion"][subname]["uptake_rate"] = sval
        #     self.celldef_tab.param_d[cdname]["secretion"][subname]["net_export_rate"] = sval

        self.current_substrate = subname
        # self.substrate_name.setText(subname)

        # item_idx = self.tree.indexFromItem(self.tree.currentItem()).row() 

        num_items = self.tree.invisibleRootItem().childCount()
        # print("tree has num_items = ",num_items)
        treeitem = QTreeWidgetItem([subname])
        treeitem.setFlags(treeitem.flags() | QtCore.Qt.ItemIsEditable)
        self.tree.insertTopLevelItem(num_items,treeitem)
        self.tree.setCurrentItem(treeitem)

        self.tree_item_clicked_cb(treeitem, 0)

    #----------------------------------------------------------------------
    # @QtCore.Slot()
    def copy_substrate(self):
        # print('------ copy_substrate')
        subname = "substrate%02d" % self.new_substrate_count
        # Make a new substrate (that's a copy of the currently selected one)
        # self.param_d[subname] = self.param_d[self.current_substrate].copy()  #rwh - "copy()" is critical
        self.param_d[subname] = copy.deepcopy(self.param_d[self.current_substrate])
        self.param_d[subname]["name"] = subname


        # for k in self.param_d.keys():
        #     print(" ===>>> ",k, " : ", self.param_d[k])

        self.new_substrate_count += 1

        # self.celldef_tab.add_new_substrate_comboboxes(subname)

        self.celldef_tab.add_new_substrate(subname)

        self.current_substrate = subname
        # self.substrate_name.setText(subname)

        # item_idx = self.tree.indexFromItem(self.tree.currentItem()).row() 

        num_items = self.tree.invisibleRootItem().childCount()
        # print("tree has num_items = ",num_items)
        treeitem = QTreeWidgetItem([subname])
        treeitem.setFlags(treeitem.flags() | QtCore.Qt.ItemIsEditable)
        self.tree.insertTopLevelItem(num_items,treeitem)
        self.tree.setCurrentItem(treeitem)

        self.tree_item_clicked_cb(treeitem, 0)
        
    #----------------------------------------------------------------------
    def show_delete_warning(self):
        msgBox = QMessageBox()
        msgBox.setIcon(QMessageBox.Information)
        msgBox.setText("Not allowed to delete all substrates.")
        #    msgBox.setWindowTitle("Example")
        # msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        msgBox.setStandardButtons(QMessageBox.Ok)
        # msgBox.buttonClicked.connect(msgButtonClick)

        returnValue = msgBox.exec()
        if returnValue == QMessageBox.Ok:
            print('OK clicked')

    #----------------------------------------------------------------------
    # @QtCore.Slot()
    def delete_substrate(self):
        num_items = self.tree.invisibleRootItem().childCount()
        print('------ delete_substrate: num_items=',num_items)
        if num_items == 1:
            # print("Not allowed to delete all substrates.")
            # QMessageBox.information(self, "Not allowed to delete all substrates")
            self.show_delete_warning()
            return

        # rwh: BEWARE of mutating the dict?
        del self.param_d[self.current_substrate]

        # may need to make a copy instead??
        # new_dict = {key:val for key, val in self.param_d.items() if key != 'Mani'}
        # self.param_d = new_dict


        # for k in self.param_d.keys():
        #     print(" ===>>> ",k, " : ", self.param_d[k])

        item_idx = self.tree.indexFromItem(self.tree.currentItem()).row() 
        # print('------      item_idx=',item_idx)
        self.tree.takeTopLevelItem(self.tree.indexOfTopLevelItem(self.tree.currentItem()))

        self.celldef_tab.delete_substrate(item_idx)
        # print('------      new name=',self.tree.currentItem().text(0))
        self.current_substrate = self.tree.currentItem().text(0)


    #----------------------------------------------------------------------
    # @QtCore.Slot()
    # def save_xml(self):
    #     # self.text.setText(random.choice(self.hello))
    #     pass

    # When a substrate is selected(via double-click) and renamed
    def tree_item_changed_cb(self, it,col):
        # print('--------- tree_item_changed_cb():', it, col, it.text(col) )  # col=0 always

        prev_name = self.current_substrate
        # print('prev_name= ',prev_name)
        self.current_substrate = it.text(col)
        self.param_d[self.current_substrate] = self.param_d.pop(prev_name)  # sweet
        # print('self.current_substrate= ',self.current_substrate )
        # for k in self.param_d.keys():
        #     print(" ===>>> ",k, " : ", self.param_d[k])
        # print()

        self.celldef_tab.renamed_substrate(prev_name, self.current_substrate)

    #----------------------------------------------------------------------
    def tree_item_sel_changed_cb(self, it,col):
        # print('--------- tree_item_sel_changed_cb():', it, col, it.text(col) )  # col=0 always

        prev_current_substrate = self.current_substrate
        self.current_substrate = it.text(col)
        self.param_d[self.current_substrate] = self.param_d.pop(prev_current_substrate)  # sweet
        # print('self.current_substrate= ',self.current_substrate )
        # for k in self.param_d.keys():
        #     print(" ===>>> ",k, " : ", self.param_d[k])
        # print()

    # def tree_item_changed_cb2(self, it,col):
    #     print('--------- tree_item_changed_cb2():', it, col, it.text(col) )  # col=0 always

    #     self.current_substrate = it.text(col)
    #     print('self.current_substrate= ',self.current_substrate )


    #----------------------------------------------------------------------
    # Update the widget values with values from param_d
    def tree_item_clicked_cb(self, it,col):
        # print('--------- tree_item_clicked_cb():', it, col, it.text(col) )  # col=0 always
        self.current_substrate = it.text(col)
        # print('self.current_substrate= ',self.current_substrate )
        # print('self.= ',self.tree.indexFromItem )

        # self.param.clear()

        # fill in the GUI with this one's params
        # self.fill_gui(self.current_substrate)

        # self.substrate_name.setText(self.param_d[self.current_substrate]["name"])
        self.diffusion_coef.setText(self.param_d[self.current_substrate]["diffusion_coef"])
        self.decay_rate.setText(self.param_d[self.current_substrate]["decay_rate"])
        self.init_cond.setText(self.param_d[self.current_substrate]["init_cond"])
        self.dirichlet_bc.setText(self.param_d[self.current_substrate]["dirichlet_bc"])
        self.dirichlet_bc_enabled.setChecked(self.param_d[self.current_substrate]["dirichlet_enabled"])

        # xmin = self.param_d[self.current_substrate]["dirichlet_xmin"]
        # print("    xmin=",xmin)
        if self.dirichlet_options_exist:
            val = self.param_d[self.current_substrate]["dirichlet_xmin"]
            # print('--------- tree_item_clicked_cb(): dirichlet_xmin=', val)
            self.dirichlet_xmin.setText(val)
            self.dirichlet_xmax.setText(self.param_d[self.current_substrate]["dirichlet_xmax"])
            self.dirichlet_ymin.setText(self.param_d[self.current_substrate]["dirichlet_ymin"])
            self.dirichlet_ymax.setText(self.param_d[self.current_substrate]["dirichlet_ymax"])
            self.dirichlet_zmin.setText(self.param_d[self.current_substrate]["dirichlet_zmin"])
            self.dirichlet_zmax.setText(self.param_d[self.current_substrate]["dirichlet_zmax"])

            # QCheckBoxs
            self.enable_xmin.setChecked(self.param_d[self.current_substrate]["enable_xmin"])
            self.enable_xmax.setChecked(self.param_d[self.current_substrate]["enable_xmax"])
            self.enable_ymin.setChecked(self.param_d[self.current_substrate]["enable_ymin"])
            self.enable_ymax.setChecked(self.param_d[self.current_substrate]["enable_ymax"])
            self.enable_zmin.setChecked(self.param_d[self.current_substrate]["enable_zmin"])
            self.enable_zmax.setChecked(self.param_d[self.current_substrate]["enable_zmax"])


        # global to all substrates
        self.gradients.setChecked(self.param_d["gradients"])
        self.track_in_agents.setChecked(self.param_d["track_in_agents"])


    #----------------------------------------------------------------------
# 		<variable name="substrate" units="dimensionless" ID="0">
# 			<physical_parameter_set>
# 				<diffusion_coefficient units="micron^2/min">100000.0</diffusion_coefficient>
# 				<decay_rate units="1/min">10</decay_rate>  
# 			</physical_parameter_set>
# 			<initial_condition units="mmHg">0</initial_condition>
# 			<Dirichlet_boundary_condition units="mmHg" enabled="true">0</Dirichlet_boundary_condition>
# <!-- use this block to set Dirichlet boundary conditions on individual boundaries --> 
# <!--
# 			<Dirichlet_options>
# 				<boundary_value ID="xmin" enabled="false">0</boundary_value>
# 				<boundary_value ID="xmax" enabled="false">0</boundary_value>
# 				<boundary_value ID="ymin" enabled="false">0</boundary_value>
# 				<boundary_value ID="ymax" enabled="false">0</boundary_value>
# 				<boundary_value ID="zmin" enabled="false">1</boundary_value>
# 				<boundary_value ID="zmax" enabled="false">0</boundary_value>
# 			</Dirichlet_options>
# -->
#  		</variable>
    def populate_tree(self):
        print("=======================  microenv populate_tree  ======================= ")
        uep = self.xml_root.find(".//microenvironment_setup")
        if uep:
            # self.substrate.clear()
            # self.param[substrate_name] = {}  # a dict of dicts

            self.tree.clear()
            idx = 0
            # <microenvironment_setup>
		    #   <variable name="food" units="dimensionless" ID="0">
            for var in uep:
                # print(cell_def.attrib['name'])
                if var.tag == 'variable':
                    substrate_name = var.attrib['name']
                    self.current_substrate = substrate_name  # do this for the callback methods (rf. BEWARE below)
                    if idx == 0:
                        # self.current_substrate = substrate_name
                        substrate_0th = substrate_name
                    self.param_d[substrate_name] = {}

                    # self.param_d[substrate_name]["name"] = substrate_name

                    treeitem = QTreeWidgetItem([substrate_name])
                    treeitem.setFlags(treeitem.flags() | QtCore.Qt.ItemIsEditable)

                    # self.substrate[var_name] = {}  # a dict of dicts
                    self.tree.insertTopLevelItem(idx,treeitem)
                    if idx == 0:  # select the 1st (0th) entry
                        self.tree.setCurrentItem(treeitem)

                    # Now fill the param dict for each substrate and the Qt widget values for the 0th

                    idx += 1

                    var_param_path = self.xml_root.find(".//microenvironment_setup//variable[" + str(idx) + "]//physical_parameter_set")
                    var_path = self.xml_root.find(".//microenvironment_setup//variable[" + str(idx) + "]")

                    # self.substrate_name.setText(var.attrib['name'])
                    diffusion_coef = var_param_path.find('.//diffusion_coefficient').text
                    # self.substrate["diffusion_coef"] = diffusion_coef
                    self.param_d[substrate_name]["diffusion_coef"] = diffusion_coef
                    if idx == 1:
                        self.diffusion_coef.setText(diffusion_coef)

                    decay_rate = var_param_path.find('.//decay_rate').text
                    # self.substrate["decay_rate"] = decay_rate
                    self.param_d[substrate_name]["decay_rate"] = decay_rate
                    if idx == 1:
                        self.decay_rate.setText(decay_rate)

                    init_cond = var_path.find('.//initial_condition').text
                    # self.substrate["init_cond"] = init_cond
                    self.param_d[substrate_name]["init_cond"] = init_cond
                    if idx == 1:
                        self.init_cond.setText(init_cond)

			# <Dirichlet_boundary_condition units="dimensionless" enabled="false">1</Dirichlet_boundary_condition>
                    dirichlet_bc_path = var_path.find('.//Dirichlet_boundary_condition')
                    dirichlet_bc = dirichlet_bc_path.text
                    # self.substrate["init_cond"] = init_cond
                    self.param_d[substrate_name]["dirichlet_bc"] = dirichlet_bc
                    # if idx == 1:
                    #     self.dirichlet_bc.setText(dirichlet_bc)

                    if "false" in dirichlet_bc_path.attrib['enabled'].lower():
                        self.param_d[substrate_name]["dirichlet_enabled"] = False
                        # if idx == 1:
                        #     self.dirichlet_bc_enabled.setChecked(False)
                    else:
                        self.param_d[substrate_name]["dirichlet_enabled"] = True
                        # if idx == 1:
                        #     self.dirichlet_bc_enabled.setChecked(True)
                        # self.dirichlet_bc_enabled.setChecked(self.param_d[self.current_substrate]["dirichlet_enabled"])

                    # 			<Dirichlet_options>
                    # 				<boundary_value ID="xmin" enabled="false">0</boundary_value>
                    # 				<boundary_value ID="xmax" enabled="false">0</boundary_value>

                    self.param_d[substrate_name]["dirichlet_xmin"] = "0"
                    self.param_d[substrate_name]["dirichlet_xmax"] = "0"
                    self.param_d[substrate_name]["dirichlet_ymin"] = "0"
                    self.param_d[substrate_name]["dirichlet_ymax"] = "0"
                    self.param_d[substrate_name]["dirichlet_zmin"] = "0"
                    self.param_d[substrate_name]["dirichlet_zmax"] = "0"
                    self.param_d[substrate_name]["enable_xmin"] = False
                    self.param_d[substrate_name]["enable_xmax"] = False
                    self.param_d[substrate_name]["enable_ymin"] = False
                    self.param_d[substrate_name]["enable_ymax"] = False
                    self.param_d[substrate_name]["enable_zmin"] = False
                    self.param_d[substrate_name]["enable_zmax"] = False

                    self.dirichlet_options_exist = True  # rwh/todo - how to handle this?
                    options_path = var_path.find('.//Dirichlet_options')
                    if options_path:
                        # self.dirichlet_options_exist = True
                        for bv in options_path:
                            print("bv = ",bv)
                            if "xmin" in bv.attrib['ID'].lower():
                                self.param_d[substrate_name]["dirichlet_xmin"] = bv.text
                                print("   -------- ",substrate_name, ":  dirichlet_xmin = ",bv.text)

                                # BEWARE: doing a 'setText' here will invoke the callback associated with
                                # the widget (e.g., self.dirichlet_xmin.textChanged.connect(self.dirichlet_xmin_changed))
                                # if idx == 1:
                                #     self.dirichlet_xmin.setText(bv.text)

                                if "true" in bv.attrib['enabled'].lower():
                                    # self.param_d[self.current_substrate]["enable_xmin"] = True
                                    self.param_d[substrate_name]["enable_xmin"] = True
                            elif "xmax" in bv.attrib['ID']:
                                self.param_d[substrate_name]["dirichlet_xmax"] = bv.text
                                # if idx == 1:
                                #     self.dirichlet_xmax.setText(bv.text)
                                if "true" in bv.attrib['enabled'].lower():
                                    self.param_d[substrate_name]["enable_xmax"] = True
                            elif "ymin" in bv.attrib['ID']:
                                self.param_d[substrate_name]["dirichlet_ymin"] = bv.text
                                # if idx == 1:
                                #     self.dirichlet_ymin.setText(bv.text)
                                if "true" in bv.attrib['enabled'].lower():
                                    self.param_d[substrate_name]["enable_ymin"] = True
                            elif "ymax" in bv.attrib['ID']:
                                self.param_d[substrate_name]["dirichlet_ymax"] = bv.text
                                self.dirichlet_ymax.setText(bv.text)
                                if "true" in bv.attrib['enabled'].lower():
                                    self.param_d[substrate_name]["enable_ymax"] = True
                            elif "zmin" in bv.attrib['ID']:
                                self.param_d[substrate_name]["dirichlet_zmin"] = bv.text
                                # self.dirichlet_zmin.setText(bv.text)
                                if "true" in bv.attrib['enabled'].lower():
                                    self.param_d[substrate_name]["enable_zmin"] = True
                            elif "zmax" in bv.attrib['ID']:
                                self.param_d[substrate_name]["dirichlet_zmax"] = bv.text
                                # self.dirichlet_zmax.setText(bv.text)
                                if "true" in bv.attrib['enabled'].lower():
                                    self.param_d[substrate_name]["enable_zmax"] = True
                    else:
                        # self.dirichlet_options_exist = False
                        self.param_d[substrate_name]["enable_xmin"] = False
                        self.param_d[substrate_name]["enable_xmax"] = False
                        self.param_d[substrate_name]["enable_ymin"] = False
                        self.param_d[substrate_name]["enable_ymax"] = False
                        self.param_d[substrate_name]["enable_zmin"] = False
                        self.param_d[substrate_name]["enable_zmax"] = False

            # </variable>
            # <options>
            # 	<calculate_gradients>true</calculate_gradients>
            # 	<track_internalized_substrates_in_each_agent>false</track_internalized_substrates_in_each_agent>
                elif var.tag == 'options':
                    self.param_d["gradients"] = False
                    self.param_d["track_in_agents"] = False
                    # self.gradients.setChecked(False)
                    # self.track_in_agents.setChecked(False)
                    for opt in var:
                        print("------- options: ",opt)
                        if "calculate_gradients" in opt.tag:
                            if "true" in opt.text.lower():
                                # self.gradients.setChecked(True)
                                self.param_d["gradients"] = True
                        elif "track_internalized_substrates_in_each_agent" in opt.tag:
                            if "true" in opt.text.lower():
                                # self.track_in_agents.setChecked(True)
                                self.param_d["track_in_agents"] = True

                    

            # options_path = uep.find(".//options")
            # print(" ---- options_path = ", options_path)
            # gradients_path = options_path.find(".//calculate_gradients")
            # # gradients_path = options_path.find("calculate_gradients")
            # print(" ---- gradients_path = ", gradients_path)
            # print(" ---- gradients_path.tag = ", gradients_path.tag)
            # print(" ---- gradients_path.text = ", gradients_path.text)
            # if "true" in gradients_path.text.lower():
            #     print(" found: gradients_path ...//calculate_gradients = true")
            #     self.param_d[self.current_substrate]["gradients"] = True

            # track_path = options_path.find(".//track_internalized_substrates_in_each_agent")
            # print(" ---- track_path.text = ", track_path.text)
            # print(" ---- track_path = ", track_path)
            # if track_path:
            #     print(" found: track_path ...//track_internalized_substrates_in_each_agent")
            # if "true" in track_path.text.lower():
            #     print(" found: track_path  = true")
            #     self.param_d[self.current_substrate]["track_in_agents"] = True

        self.current_substrate = substrate_0th
        self.tree.setCurrentItem(self.tree.topLevelItem(0))  # select the top (0th) item
        self.tree_item_clicked_cb(self.tree.topLevelItem(0), 0)  # and invoke its callback to fill widget values

        print("\n\n=======================  leaving microenv populate_tree  ======================= ")
        # for k in self.param_d.keys():
        #     print(" ===>>> ",k, " : ", self.param_d[k])

#        ---- populate_tree(): self.param_d =  {'director signal': {'diffusion_coef': '1000', 'decay_rate': '.4', 'init_cond': '0', 'dirichlet_bc': '1', 'dirichlet_enabled': False, 'enable_xmin': False, 'enable_xmax': False, 'enable_ymin': False, 'enable_ymax': False, 'enable_zmin': False, 'enable_zmax': False, 'dirichlet_xmin': '-11', 'dirichlet_xmax': '11', 'dirichlet_ymin': '-12', 'dirichlet_ymax': '12', 'dirichlet_zmin': '-13', 'dirichlet_zmax': '13'}, 'cargo signal': {'diffusion_coef': '1000', 'decay_rate': '.4', 'init_cond': '0', 'dirichlet_bc': '1', 'dirichlet_enabled': False, 'enable_xmin': False, 'enable_xmax': False, 'enable_ymin': False, 'enable_ymax': False, 'enable_zmin': False, 'enable_zmax': False, 'dirichlet_xmin': '-11', 'dirichlet_xmax': '11', 'dirichlet_ymin': '-12', 'dirichlet_ymax': '12', 'dirichlet_zmin': '-13', 'dirichlet_zmax': '13'}}


    #----------------------------------------------------------------------------
    def first_substrate_name(self):
        uep = self.xml_root.find(".//microenvironment_setup//variable")
        if uep:
                return(uep.attrib['name'])


            #----------------------------------------------------------------------------
            # Read values from the params_d and generate XML

            # 	<microenvironment_setup>
            # 	<variable name="director signal" units="dimensionless" ID="0">
            # 		<physical_parameter_set>
            # 			<diffusion_coefficient units="micron^2/min">1000</diffusion_coefficient>
            # 			<decay_rate units="1/min">.1</decay_rate>  
            # 		</physical_parameter_set>
            # 		<initial_condition units="dimensionless">0</initial_condition>
            # 		<Dirichlet_boundary_condition units="dimensionless" enabled="false">1</Dirichlet_boundary_condition>
            # 	</variable>
                
            # 	<variable name="cargo signal" units="dimensionless" ID="1">
            # 		<physical_parameter_set>
            # 			<diffusion_coefficient units="micron^2/min">1000</diffusion_coefficient>
            # 			<decay_rate units="1/min">.4</decay_rate>  
            # 		</physical_parameter_set>
            # 		<initial_condition units="dimensionless">0</initial_condition>
            # 		<Dirichlet_boundary_condition units="dimensionless" enabled="false">1</Dirichlet_boundary_condition>
            # 	</variable>
                
            # 	<options>
            # 		<calculate_gradients>true</calculate_gradients>
            # 		<track_internalized_substrates_in_each_agent>false</track_internalized_substrates_in_each_agent>
                    
            # 		<initial_condition type="matlab" enabled="false">
            # 			<filename>./config/initial.mat</filename>
            # 		</initial_condition>
                    
            # 		<dirichlet_nodes type="matlab" enabled="false">
            # 			<filename>./config/dirichlet.mat</filename>
            # 		</dirichlet_nodes>
            # 	</options>
            # </microenvironment_setup>


    # <microenvironment_setup>
	# 	<variable name="oxygen" units="mmHg" ID="0">
	# 		<physical_parameter_set>
	# 			<diffusion_coefficient units="micron^2/min">421.0</diffusion_coefficient>
	# 			<decay_rate units="1/min">.41</decay_rate>  
	# 		</physical_parameter_set>
	# 		<initial_condition units="mmHg">41.0</initial_condition>
	# 		<Dirichlet_boundary_condition units="mmHg" enabled="true">41.1</Dirichlet_boundary_condition>
    #         <Dirichlet_options>
 	# 			<boundary_value ID="xmin" enabled="false">1</boundary_value>
 	# 			<boundary_value ID="xmax" enabled="true">2</boundary_value>
 	# 			<boundary_value ID="ymin" enabled="false">3</boundary_value>
 	# 			<boundary_value ID="ymax" enabled="true">4</boundary_value>
 	# 			<boundary_value ID="zmin" enabled="false">5</boundary_value>
 	# 			<boundary_value ID="zmax" enabled="true">6</boundary_value>
 	# 		</Dirichlet_options>
	# 	</variable>
	
	# 	<variable name="glue" units="dimensionless" ID="1">
	# 		<physical_parameter_set>
	# 			<diffusion_coefficient units="micron^2/min">422.0</diffusion_coefficient>
	# 			<decay_rate units="1/min">.42</decay_rate>  
	# 		</physical_parameter_set>
	# 		<initial_condition units="mmHg">42.0</initial_condition>
	# 		<Dirichlet_boundary_condition units="mmHg" enabled="false">42.1</Dirichlet_boundary_condition>
	# 	</variable>
		
	# 	<options>
	# 		<calculate_gradients>true</calculate_gradients>
	# 		<track_internalized_substrates_in_each_agent>false</track_internalized_substrates_in_each_agent>
			 
	# 		<initial_condition type="matlab" enabled="false">
	# 			<filename>./config/initial.mat</filename>
	# 		</initial_condition>
			 
	# 		<dirichlet_nodes type="matlab" enabled="false">
	# 			<filename>./config/dirichlet.mat</filename>
	# 		</dirichlet_nodes>
	# 	</options>
	# </microenvironment_setup>	

    def iterate_tree(self, node, count, subs):
        for idx in range(count):
            item = node.child(idx)
            # print('******* State: %s, Text: "%s"' % (Item.checkState(3), Item.text(0)))
            subs.append(item.text(0))
            child_count = item.childCount()
            if child_count > 0:
                self.iterate_tree(item, child_count)
                
    def fill_xml(self):
        print("----------- microenv_tab.py: fill_xml(): ----------")
        uep = self.xml_root.find('.//microenvironment_setup') # guaranteed to exist since we start with a valid model
        vp = []   # pointers to <variable> nodes
        if uep:
            # Begin by removing all previously defined substrates in the .xml
            for var in uep.findall('variable'):
                uep.remove(var)
                # vp.append(var)
        # self.tree_status()

        # Obtain a list of all substrates in self.tree (QTreeWidget()). Used below.
        substrates_in_tree = []
        num_subs = self.tree.invisibleRootItem().childCount()  # rwh: get number of items in tree
        print('num subtrates = ',num_subs)
        self.iterate_tree(self.tree.invisibleRootItem(), num_subs, substrates_in_tree)
        print("substrates_in_tree =",substrates_in_tree)

        uep = self.xml_root.find('.//microenvironment_setup')
        indent1 = '\n'
        indent6 = '\n      '
        indent8 = '\n        '
        indent10 = '\n          '

        idx = 0
        for substrate in self.param_d.keys():
            print('key in param_d.keys() = ',substrate)
            if substrate in substrates_in_tree:
                print("matched! ",substrate)
	# 	<variable name="glue" units="dimensionless" ID="1">
	# 		<physical_parameter_set>
	# 			<diffusion_coefficient units="micron^2/min">422.0</diffusion_coefficient>
	# 			<decay_rate units="1/min">.42</decay_rate>  
	# 		</physical_parameter_set>
	# 		<initial_condition units="mmHg">42.0</initial_condition>
	# 		<Dirichlet_boundary_condition units="mmHg" enabled="false">42.1</Dirichlet_boundary_condition>
                # elm = ET.Element(substrate)
                # elm = ET.Element(substrate+'\n', {'foo':'bar'})


        # self.param_d[self.current_substrate]["diffusion_coef"] = text
        # self.param_d[self.current_substrate]["decay_rate"] = text
        # self.param_d[self.current_substrate]["init_cond"] = text
        # self.param_d[self.current_substrate]["dirichlet_bc"] = text
        # self.param_d[self.current_substrate]["dirichlet_enabled"] = self.dirichlet_bc_enabled.isChecked()
                elm = ET.Element("variable", 
                        {"name":substrate, "units":"dimensionless", "ID":str(idx)})
                elm.tail = '\n' + indent6
                elm.text = indent8
                subelm = ET.SubElement(elm, 'physical_parameter_set')
                subelm.text = indent10
                subelm.tail = indent8
                subelm2 = ET.SubElement(subelm, "diffusion_coefficient",{"units":"micron^2/min"})
                subelm2.text = self.param_d[substrate]["diffusion_coef"]
                subelm2.tail = indent10
                subelm2 = ET.SubElement(subelm, "decay_rate",{"units":"1/min"})
                subelm2.text = self.param_d[substrate]["decay_rate"]
                subelm2.tail = indent8

                subelm = ET.SubElement(elm, 'initial_condition', {"units":"mmHg"})
                subelm.text = self.param_d[substrate]["init_cond"]
                subelm.tail = indent8
                subelm = ET.SubElement(elm, "Dirichlet_boundary_condition",
                        {"units":"mmHg", "enabled":str(self.param_d[self.current_substrate]["dirichlet_enabled"])})
                subelm.text = self.param_d[substrate]["dirichlet_bc"]
                subelm.tail = indent6
                        
                #              {'text':"foo",
                #               'xmlUrl':"bar",
                #               'htmlUrl':"grrr",
                #               })
                # uep.append(elm)
                uep.insert(idx,elm)
                idx += 1

        # print(prettify(self.xml_root))

	# 	<variable name="oxygen" units="mmHg" ID="0">
	# 		<physical_parameter_set>
	# 			<diffusion_coefficient units="micron^2/min">421.0</diffusion_coefficient>
	# 			<decay_rate units="1/min">.41</decay_rate>  
	# 		</physical_parameter_set>
	# 		<initial_condition units="mmHg">41.0</initial_condition>
	# 		<Dirichlet_boundary_condition units="mmHg" enabled="true">41.1</Dirichlet_boundary_condition>
    #         <Dirichlet_options>
 	# 			<boundary_value ID="xmin" enabled="false">1</boundary_value>
 	# 			<boundary_value ID="xmax" enabled="true">2</boundary_value>
 	# 			<boundary_value ID="ymin" enabled="false">3</boundary_value>
 	# 			<boundary_value ID="ymax" enabled="true">4</boundary_value>
 	# 			<boundary_value ID="zmin" enabled="false">5</boundary_value>
 	# 			<boundary_value ID="zmax" enabled="true">6</boundary_value>
 	# 		</Dirichlet_options>
	# 	</variable>


        # ------ Finally, append the flags that apply to all substrates
        if self.gradients.isChecked():
            self.xml_root.find(".//options//calculate_gradients").text = 'true'
        else:
            self.xml_root.find(".//options//calculate_gradients").text = 'false'

        if self.track_in_agents.isChecked():
            self.xml_root.find(".//options//track_internalized_substrates_in_each_agent").text = 'true'
        else:
            self.xml_root.find(".//options//track_internalized_substrates_in_each_agent").text = 'false'
    
    def clear_gui(self):
        pass
    def create_main_tab(self):
        main_tab = QWidget()
        glayout = QGridLayout()
        #------------------
        self.microenv_tabs_layout = QGridLayout()
        self.microenv_tabs_layout.addWidget(self.tab_widget, 0,0,1,1) # w, row, column, rowspan, colspan
        #------------------
        hbox = QHBoxLayout()
        label = QLabel("diffusion coefficient")
        label.setFixedWidth(self.label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.diffusion_coef = QLineEdit()
        self.diffusion_coef.setValidator(QtGui.QDoubleValidator())
        self.diffusion_coef.textChanged.connect(self.diffusion_coef_changed)
        # self.diffusion_coef.enter.connect(self.save_xml)
        hbox.addWidget(self.diffusion_coef)

        units = QLabel("micron^2/min")
        units.setFixedWidth(self.units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)

        #----------
        hbox = QHBoxLayout()
        label = QLabel("decay rate")
        label.setFixedWidth(self.label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.decay_rate = QLineEdit()
        self.decay_rate.setValidator(QtGui.QDoubleValidator())
        self.decay_rate.textChanged.connect(self.decay_rate_changed)
        # self.decay_rate.enter.connect(self.save_xml)
        hbox.addWidget(self.decay_rate)

        units = QLabel("1/min")
        units.setFixedWidth(self.units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)

        #----------
        hbox = QHBoxLayout()
        label = QLabel("initial condition")
        label.setFixedWidth(self.label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.init_cond = QLineEdit()
        self.init_cond.setValidator(QtGui.QDoubleValidator())
        self.init_cond.textChanged.connect(self.init_cond_changed)
        # self.init_cond.enter.connect(self.save_xml)
        hbox.addWidget(self.init_cond)

        units = QLabel("mmol")
        units.setFixedWidth(self.units_width)
        hbox.addWidget(units)
        self.vbox.addLayout(hbox)
        #----------

        hbox = QHBoxLayout()
        label = QLabel("Dirichlet BC")
        label.setFixedWidth(self.label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.dirichlet_bc = QLineEdit()
        self.dirichlet_bc.setValidator(QtGui.QDoubleValidator())
        self.dirichlet_bc.textChanged.connect(self.dirichlet_bc_changed)
        # self.bdy_cond.enter.connect(self.save_xml)
        hbox.addWidget(self.dirichlet_bc)

        units = QLabel("mmol")
        units.setFixedWidth(self.units_width)
        hbox.addWidget(units)

        self.dirichlet_bc_enabled = QCheckBox("on")
        self.dirichlet_bc_enabled.stateChanged.connect(self.dirichlet_toggle_cb)
        # self.motility_enabled.setAlignment(QtCore.Qt.AlignRight)
        # label.setFixedWidth(self.label_width)
        hbox.addWidget(self.dirichlet_bc_enabled)

        self.vbox.addLayout(hbox)

        #--------------------------

        dirichlet_options_bdy = QLabel("Dirichlet options per boundary:")
        # units.setFixedWidth(self.units_width)
        self.vbox.addWidget(dirichlet_options_bdy)

        #----
        hbox = QHBoxLayout()
        label = QLabel("xmin:")
        label.setFixedWidth(self.label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.dirichlet_xmin = QLineEdit()
        self.dirichlet_xmin.setValidator(QtGui.QDoubleValidator())
        self.dirichlet_xmin.textChanged.connect(self.dirichlet_xmin_changed)
        hbox.addWidget(self.dirichlet_xmin)

        self.enable_xmin = QCheckBox("on")
        self.enable_xmin.stateChanged.connect(self.enable_xmin_cb)
        # self.motility_enabled.setAlignment(QtCore.Qt.AlignRight)
        # label.setFixedWidth(self.label_width)
        hbox.addWidget(self.enable_xmin)
        self.vbox.addLayout(hbox)
        #----
        hbox = QHBoxLayout()
        label = QLabel("xmax:")
        label.setFixedWidth(self.label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.dirichlet_xmax = QLineEdit()
        self.dirichlet_xmax.setValidator(QtGui.QDoubleValidator())
        self.dirichlet_xmax.textChanged.connect(self.dirichlet_xmax_changed)
        hbox.addWidget(self.dirichlet_xmax)

        self.enable_xmax = QCheckBox("on")
        self.enable_xmax.stateChanged.connect(self.enable_xmax_cb)
        hbox.addWidget(self.enable_xmax)
        self.vbox.addLayout(hbox)
        #---------
        hbox = QHBoxLayout()
        label = QLabel("ymin:")
        label.setFixedWidth(self.label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.dirichlet_ymin = QLineEdit()
        self.dirichlet_ymin.setValidator(QtGui.QDoubleValidator())
        self.dirichlet_ymin.textChanged.connect(self.dirichlet_ymin_changed)
        hbox.addWidget(self.dirichlet_ymin)

        self.enable_ymin = QCheckBox("on")
        self.enable_ymin.stateChanged.connect(self.enable_ymin_cb)
        # self.motility_enabled.setAlignment(QtCore.Qt.AlignRight)
        # label.setFixedWidth(self.label_width)
        hbox.addWidget(self.enable_ymin)
        self.vbox.addLayout(hbox)
        #----
        hbox = QHBoxLayout()
        label = QLabel("ymax:")
        label.setFixedWidth(self.label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.dirichlet_ymax = QLineEdit()
        self.dirichlet_ymax.setValidator(QtGui.QDoubleValidator())
        self.dirichlet_ymax.textChanged.connect(self.dirichlet_ymax_changed)
        hbox.addWidget(self.dirichlet_ymax)

        self.enable_ymax = QCheckBox("on")
        self.enable_ymax.stateChanged.connect(self.enable_ymax_cb)
        hbox.addWidget(self.enable_ymax)
        self.vbox.addLayout(hbox)
        #---------
        hbox = QHBoxLayout()
        label = QLabel("zmin:")
        label.setFixedWidth(self.label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.dirichlet_zmin = QLineEdit()
        self.dirichlet_zmin.setValidator(QtGui.QDoubleValidator())
        self.dirichlet_zmin.textChanged.connect(self.dirichlet_zmin_changed)
        hbox.addWidget(self.dirichlet_zmin)

        self.enable_zmin = QCheckBox("on")
        self.enable_zmin.stateChanged.connect(self.enable_zmin_cb)
        # self.motility_enabled.setAlignment(QtCore.Qt.AlignRight)
        # label.setFixedWidth(self.label_width)
        hbox.addWidget(self.enable_zmin)
        self.vbox.addLayout(hbox)
        #----
        hbox = QHBoxLayout()
        label = QLabel("zmax:")
        label.setFixedWidth(self.label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        hbox.addWidget(label)

        self.dirichlet_zmax = QLineEdit()
        self.dirichlet_zmax.setValidator(QtGui.QDoubleValidator())
        self.dirichlet_zmax.textChanged.connect(self.dirichlet_zmax_changed)
        hbox.addWidget(self.dirichlet_zmax)

        self.enable_zmax = QCheckBox("on")
        self.enable_zmax.stateChanged.connect(self.enable_zmax_cb)
        hbox.addWidget(self.enable_zmax)
        self.vbox.addLayout(hbox)

        #-------------
        # Toggles for overall microenv (all substrates)
        self.vbox.addWidget(QHLine())

        hbox = QHBoxLayout()
        hbox.addWidget(QLabel("For all substrates: "))

        self.gradients = QCheckBox("calculate gradients")
        self.gradients.stateChanged.connect(self.gradients_cb)
        hbox.addWidget(self.gradients)

        self.track_in_agents = QCheckBox("track in agents")
        self.track_in_agents.stateChanged.connect(self.track_in_agents_cb)
        hbox.addWidget(self.track_in_agents)
        self.vbox.addLayout(hbox)

        self.vbox.addStretch()

        self.microenv_params.setLayout(self.vbox)

        self.scroll_area.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll_area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setWidget(self.microenv_params)

        self.layout = QVBoxLayout(self)
        main_tab.setLayout(glayout)
        return main_tab
    def create_secretion_tab(self):
        secretion_tab = QWidget()
        glayout = QGridLayout()

        label = QLabel("Phenotype: secretion")
        label.setStyleSheet("background-color: orange")
        label.setAlignment(QtCore.Qt.AlignCenter)

        self.secretion_substrate_dropdown = QComboBox()
        idr = 0
        glayout.addWidget(self.secretion_substrate_dropdown, idr,0, 1,1) # w, row, column, rowspan, colspan
        #self.secretion_substrate_dropdown.currentIndexChanged.connect(self.secretion_substrate_changed_cb)  # beware: will be triggered on a ".clear" too

        label = QLabel("secretion rate")
        label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # # label.setStyleSheet("border: 1px solid black;")
        # idr += 1
        # glayout.addWidget(label, idr,0, 1,1) # w, row, column, rowspan, colspan

        # self.secretion_rate = QLineEdit()
        # self.secretion_rate.textChanged.connect(self.secretion_rate_changed)
        # self.secretion_rate.setValidator(QtGui.QDoubleValidator())
        # glayout.addWidget(self.secretion_rate, idr,1, 1,1) # w, row, column, rowspan, colspan

        # units = QLabel("1/min")
        # units.setFixedWidth(self.units_width)
        # units.setAlignment(QtCore.Qt.AlignLeft)
        # # units.setStyleSheet("border: 1px solid black;")
        # glayout.addWidget(units, idr,2, 1,1) # w, row, column, rowspan, colspan

        # #---
        # label = QLabel("target")
        # label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # # label.setStyleSheet("border: 1px solid black;")
        # idr += 1
        # glayout.addWidget(label, idr,0, 1,1) # w, row, column, rowspan, colspan

        # self.secretion_target = QLineEdit()
        # self.secretion_target.textChanged.connect(self.secretion_target_changed)
        # self.secretion_target.setValidator(QtGui.QDoubleValidator())
        # glayout.addWidget(self.secretion_target, idr,1, 1,1) # w, row, column, rowspan, colspan

        # # units = QLabel("substrate density")
        # units = QLabel("sub. density")
        # # units.setFixedWidth(self.units_width+5)
        # # units.setFixedWidth(110)
        # units.setAlignment(QtCore.Qt.AlignLeft)
        # # units.setStyleSheet("border: 1px solid black;")
        # glayout.addWidget(units, idr,2, 1,1) # w, row, column, rowspan, colspan

        # #---
        # label = QLabel("uptake rate")
        # label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # idr += 1
        # glayout.addWidget(label, idr,0, 1,1) # w, row, column, rowspan, colspan

        # self.uptake_rate = QLineEdit()
        # self.uptake_rate.textChanged.connect(self.uptake_rate_changed)
        # self.uptake_rate.setValidator(QtGui.QDoubleValidator())
        # glayout.addWidget(self.uptake_rate, idr,1, 1,1) # w, row, column, rowspan, colspan

        # units = QLabel("1/min")
        # units.setFixedWidth(self.units_width)
        # units.setAlignment(QtCore.Qt.AlignLeft)
        # glayout.addWidget(units, idr,2, 1,1) # w, row, column, rowspan, colspan

        # #---
        # label = QLabel("net export rate")
        # label.setFixedWidth(self.label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # idr += 1
        # glayout.addWidget(label, idr,0, 1,1) # w, row, column, rowspan, colspan

        # self.secretion_net_export_rate = QLineEdit()
        # self.secretion_net_export_rate.textChanged.connect(self.secretion_net_export_rate_changed)
        # self.secretion_net_export_rate.setValidator(QtGui.QDoubleValidator())
        # glayout.addWidget(self.secretion_net_export_rate, idr,1, 1,1) # w, row, column, rowspan, colspan

        # units = QLabel("total/min")
        # units.setFixedWidth(self.units_width)
        # units.setAlignment(QtCore.Qt.AlignLeft)
        # glayout.addWidget(units, idr,2, 1,1) # w, row, column, rowspan, colspan

        # #------
        # for idx in range(11):  # rwh: hack solution to align rows
        #     blank_line = QLabel("")
        #     idr += 1
        #     glayout.addWidget(blank_line, idr,0, 1,1) # w, row, column, rowspan, colspan

        # #------
        # vlayout.setVerticalSpacing(10)  # rwh - argh
        secretion_tab.setLayout(glayout)
        return secretion_tab