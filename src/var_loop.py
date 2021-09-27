import xml.etree.ElementTree as ET

import xml.etree.ElementTree as ET
tree = ET.parse('../config/mymodel.xml')
root = tree.getroot()


up = root.find('user_parameters')

dc = up.find('PKPD_D1_central_increase_on_dose')

dc.text = str(7*float(dc.text));

tree.write('../config/mymodel_new.xml')
