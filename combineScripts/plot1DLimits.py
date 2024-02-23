import os,sys
import numpy as np
import ROOT

drawObserved = True
drawPoints = True
maskSMResonances = True

model = sys.argv[1]
limdir = sys.argv[2]
ctau = sys.argv[3]

massl = []
obsl  = []
expl  = []
m2sl  = []
m1sl  = []
p1sl  = []
p2sl  = []

fin = open("%s/limits_%s.txt"%(limdir,model),"r")
for l in fin.readlines():
    if l.startswith("#"):
        continue
    ls = l.split(",")
    if ls[0]!=model:
        continue
    if int(ls[2])!=int(ctau):
        continue
    massl.append(float(ls[1]))
    obsl .append(float(ls[3]))
    expl .append(float(ls[4]))
    m2sl .append(float(ls[5]))
    m1sl .append(float(ls[6]))
    p1sl .append(float(ls[7]))
    p2sl .append(float(ls[8]))
fin.close()

massv = np.array(massl,"d")
obsv  = np.array(obsl ,"d")
expv  = np.array(expl ,"d")
m1sv  = np.array(m1sl ,"d")
p1sv  = np.array(p1sl ,"d")

print(obsl)

miny = 999.
maxy = -1.0
if min(obsl)<miny:
    miny=min(obsl)*0.1
if min(m1sl)<miny:
    miny=min(m1sl)*0.1

if max(obsl)>maxy:
    maxy=max(obsl)*10.0
if max(p1sl)>maxy:
    maxy=max(p1sl)*10.0

#miny=0.1
#maxy=50.

gobs = ROOT.TGraph(len(massv),massv,obsv)
gobs.SetLineColor(1)
gobs.SetMarkerColor(1)
gobs.SetMarkerStyle(20)
gobs.SetLineStyle(1)
gobs.SetLineWidth(2)

gexp = ROOT.TGraph(len(massv),massv,expv)
gexp.SetLineColor(2)
gexp.SetMarkerColor(2)
gexp.SetMarkerStyle(20)
gexp.SetLineStyle(1)
gexp.SetLineWidth(2)

gm1s = ROOT.TGraph(len(massv),massv,m1sv)
gm1s.SetLineColor(2)
gm1s.SetMarkerColor(2)
gm1s.SetMarkerStyle(20)
gm1s.SetLineStyle(2)
gm1s.SetLineWidth(2)

gp1s = ROOT.TGraph(len(massv),massv,p1sv)
gp1s.SetLineColor(2)
gp1s.SetMarkerColor(2)
gp1s.SetMarkerStyle(20)
gp1s.SetLineStyle(2)
gp1s.SetLineWidth(2)

leg = ROOT.TLegend(0.6,0.65,0.87,0.85)
#leg.SetLineColor(0)
leg.SetFillColor(0)
leg.SetTextSize(0.03)
modelLeg = model
if model=="HTo2ZdTo2mu2x":
    modelLeg="h#rightarrow Z_{D}Z_{D}"
elif model=="DY3":
    modelLeg="DY_{3}"
elif model=="DYp3":
    modelLeg="DY'_{3}"
elif model=="B3mL2":
    modelLeg="B_{3}-L_{2}"
leg.SetHeader("%s (c#tau = %s mm)"%(modelLeg,ctau))
if drawObserved:
    leg.AddEntry(gobs,"Observed","L")
leg.AddEntry(gexp,"Expected","L")
leg.AddEntry(gm1s,"#pm1#sigma expected","L")

haxis = ROOT.TH1D("haxis","",10,massv[0]-0.05*massv[0],massv[len(massv)-1]+0.05*massv[len(massv)-1])
haxis.GetXaxis().SetTitle("Mass [GeV]")
haxis.GetYaxis().SetTitle("95% CL upper limit on #sigma/#sigma_{theory}")
haxis.GetXaxis().SetLabelSize(0.03)
haxis.GetYaxis().SetLabelSize(0.03)
haxis.SetMinimum(miny)
haxis.SetMaximum(maxy)
haxis.GetYaxis().SetRangeUser(miny,maxy)

line = ROOT.TLine(massv[0],1.0,massv[len(massv)-1],1.0)
line.SetLineColor(ROOT.kGray)
line.SetLineStyle(2)

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetFrameLineWidth(2)
can = ROOT.TCanvas("can","",700,600)
can.cd()
ROOT.gPad.SetTicky()
ROOT.gPad.SetLogy()
ROOT.gPad.SetLogx()

haxis.Draw()

if not drawPoints:
    gm1s.Draw("L")
    gexp.Draw("L")
    gp1s.Draw("L")
    if drawObserved:
        gobs.Draw("L")
else:
    gexp.Draw("PL")
    gm1s.Draw("L")
    gp1s.Draw("L")
    if drawObserved:
        gobs.Draw("PL")

if maskSMResonances:

    brho = ROOT.TGraph(4)
    brho.SetPoint(0, 0.7, miny+0.01*miny)
    brho.SetPoint(1, 0.89, miny+0.01*miny)
    brho.SetPoint(2, 0.89, maxy-0.05)
    brho.SetPoint(3, 0.7, maxy-0.05)
    brho.SetFillColorAlpha(ROOT.kGray, 1.0)
    brho.Draw("f")

    bphi = ROOT.TGraph(4)
    bphi.SetPoint(0, 0.91, miny+0.01*miny)
    bphi.SetPoint(1, 1.2, miny+0.01*miny)
    bphi.SetPoint(2, 1.2, maxy-0.05)
    bphi.SetPoint(3, 0.91, maxy-0.05)
    bphi.SetFillColorAlpha(ROOT.kGray, 1.0)
    bphi.Draw("f")

    bjpsi = ROOT.TGraph(4)
    bjpsi.SetPoint(0, 2.8, miny+0.01*miny)
    bjpsi.SetPoint(1, 4.1, miny+0.01*miny)
    bjpsi.SetPoint(2, 4.1, maxy-0.05)
    bjpsi.SetPoint(3, 2.8, maxy-0.05)
    bjpsi.SetFillColorAlpha(ROOT.kGray, 1.0)
    bjpsi.Draw("f")

    bups = ROOT.TGraph(4)
    bups.SetPoint(0, 8.6, miny+0.01*miny)
    bups.SetPoint(1, 11.0, miny+0.01*miny)
    bups.SetPoint(2, 11.0, maxy-0.05)
    bups.SetPoint(3, 8.6, maxy-0.05)
    bups.SetFillColorAlpha(ROOT.kGray, 1.0)
    bups.Draw("f")


line.Draw("same")

leg.Draw("same")

ROOT.gPad.RedrawAxis()
can.Update()
can.SaveAs("%s/limits_%s_ctau%s.png"%(limdir,model,ctau))
