import sys,os

CTRMAX = 1081

for n in xrange(CTRMAX+1):

    CMD = 'python drremoval.py /home/david/uboone/data/mcc84/ccpi0/larlite_reco2d_%04i.root /home/david/uboone/data/mcc84/ccpi0/proximity_%04i.root /home/david/uboone/data/mcc84/ccpi0/em_%04i.root'%(n,n,n)
    os.system(CMD)
