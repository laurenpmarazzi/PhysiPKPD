
# coding: utf-8

# In[20]:


import tellurium as te

dosingmodel_str = '''

model PKmodel()

// c1 -> central compartment
// c2 -> other peripheral compartment
// c_tumor -> tumor compartment

J0: drug -> c1; k0*drug
J1: c1 -> c2; k12*c1
J3: c2 -> c1; k21*c2
J4: c1 -> c_tumor; k34*c1
J5: c_tumor -> c1; k43*c_tumor
J6: c1 -> ; k10*c1

//Initial values

drug = 80;
c1 = 0;
c2 = 0;
c_tumor = 0;

// rate parameters

k0 = 1;
k12 = 0.0041;
k21 = 0.0174;
k34 = 0.0269;
k43 = 0.0107;
k10 =0.0394;

end
'''

rr = te.loadAntimonyModel(dosingmodel_str)
n = rr.simulate(0,600,10001,['Time','c_tumor'])
rr.plot(n)

