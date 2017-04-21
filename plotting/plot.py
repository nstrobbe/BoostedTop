#! /usr/bin/env python
import os
import glob
import math
import array
import sys
import time

import ROOT

ROOT.gROOT.ProcessLine(".L tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadLeftMargin(0.16);


def getMax(histos):
    temp_max = 0
    for h in histos:
        hmax = h.GetMaximum()
        if hmax>temp_max:
            temp_max = hmax
    return temp_max

############################################################

if __name__ == '__main__':

    filenames = [#"../out_gg_s.root",
                 #"../out_gg_a.root",
                 #"../out_pp_zp.root",
                 "../out_pp_zpbl.root",
                 "../out_pp_zpchi.root",
                 "../out_pp_zpeta.root",
                 "../out_pp_zplr.root",
                 "../out_pp_zppsi.root",
                 "../out_pp_zpssm.root"
                 ]
    labels = [f.split("/")[-1].strip(".root").lstrip("out_") for f in filenames]
    files = [ROOT.TFile(fname) for fname in filenames]
    colors = [#ROOT.kBlack,
              #ROOT.kGray+1,
              #ROOT.kRed,
              ROOT.kMagenta,
              ROOT.kBlue,
              ROOT.kGreen+2,
              ROOT.kMagenta+2,
              ROOT.kOrange,
              ROOT.kCyan+1]

    trees = [f.Get("tree") for f in files]
    
    vars = ["costheta1","costheta2","costhetastar","phi","phi1","costhetaW1","costhetaW2","phi2","phiW1","phiW2",
            "dPhi12","dPhiW1W2","dPhi1W1","dPhi2W2","mttbar"]
    vars_x = ["costheta1","costheta2","costhetastar","phi","phi1","costhetaW1","costhetaW2","phi2","phiW1","phiW2",
              "#Delta(Phi1,Phi2)","#Delta(PhiW1,PhiW2)","#Delta(Phi1,PhiW1)","#Delta(Phi2,PhiW2)","m(ttbar)"]
    vars_nbin = [20]*len(vars);
    vars_lo = [-1,  -1,   -1, -ROOT.TMath.Pi(), -ROOT.TMath.Pi(), -1,   -1, -ROOT.TMath.Pi(), -ROOT.TMath.Pi(), -ROOT.TMath.Pi(),
                0, 0, 0, 0, 0]
    vars_hi = [ 1,   1,    1,  ROOT.TMath.Pi(),  ROOT.TMath.Pi(),  1,    1,  ROOT.TMath.Pi(),  ROOT.TMath.Pi(),  ROOT.TMath.Pi(),
                ROOT.TMath.Pi(),  ROOT.TMath.Pi(),  ROOT.TMath.Pi(), ROOT.TMath.Pi(), 3000]    
    
    histos = []
    for l in labels:
        histo_l = []
        for i in range(len(vars)):
            histo_l.append( ROOT.TH1F("h_"+l+"_"+vars[i],";"+vars_x[i]+";n.d.",vars_nbin[i],vars_lo[i],vars_hi[i]) )
        histos.append(histo_l)
    
    for i_t, t in enumerate(trees):
        for i in range( t.GetEntries()):
            t.GetEntry(i)
            if i%2000 == 0: print "i = ",i
            
            for j in range(len(vars)):
                histos[i_t][j].Fill( getattr( t, vars[j] ) )

                
    ##### plot
    # scale and color
    for i_h, h in enumerate(histos):
        for j in range(len(vars)):
            h[j].Scale( 1./h[j].Integral() );
            h[j].SetLineColor( colors[i_h] )
            h[j].SetLineWidth( 3 )
        
    # on canvas
    for j in range(len(vars)):
        
        tmpcan = ROOT.TCanvas("can"+str(j),"can"+str(j),800,800);
        tmpcan.SetRightMargin(0.05)

        leg = ROOT.TLegend(0.25, 0.75, 0.92, 0.92)
        leg.SetBorderSize(0)
        leg.SetNColumns(2)
        max = getMax( [h[j] for h in histos]  )
        for i_h, h in enumerate(histos):
            h[j].SetMaximum(1.5*max)
            h[j].SetMinimum(0)
            if i_h == 0:
                h[j].Draw()
            else:
                h[j].Draw("same")
            leg.AddEntry(h[j], labels[i_h], "l")

        leg.Draw()

        tmpcan.SaveAs("plots/"+vars[j]+".eps");
        tmpcan.SaveAs("plots/"+vars[j]+".png");        
        
        
        
        
        
        
        
        
