using LightXML

xdoc = parse_file("./config/mymodel.xml")
r = root(xdoc)
up = find_element(r, "user_parameters")
d1up = find_element(up, "PKPD_D1_central_increase_on_dose")
set_content(d1up, "5e4")

mycommand = `./project ./config/mymodel.xml`
run(mycommand)