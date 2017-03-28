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


############################################################

if __name__ == '__main__':

    file_gg_s = ROOT.TFile("../out_gg_s.root");
    file_pp_zp = ROOT.TFile("../out_pp_zp.root");
    tree_gg_s = file_gg_s.Get("tree");
    tree_pp_zp = file_pp_zp.Get("tree");
    
    vars = ["costheta1","costheta2","costhetastar","phi","phi1"]
    vars_x = ["costheta1","costheta2","costhetastar","phi","phi1"]
    vars_nbin = [20]*5;
    vars_lo = [-1,   0,   -1, -ROOT.TMath.Pi(), -ROOT.TMath.Pi()]
    vars_hi = [ 1,   1,    1,  ROOT.TMath.Pi(),  ROOT.TMath.Pi()]    
    
    h_gg_s = [];
    h_pp_zp = [];
    
    for i in range(len(vars)):
        h_gg_s.append( ROOT.TH1F("h_gg_s_"+vars[i],";"+vars_x[i]+";n.d.",vars_nbin[i],vars_lo[i],vars_hi[i]) );
        h_pp_zp.append( ROOT.TH1F("h_pp_zp_"+vars[i],";"+vars_x[i]+";n.d.",vars_nbin[i],vars_lo[i],vars_hi[i]) );   

    
    for i in range( tree_gg_s.GetEntries() ):
        tree_gg_s.GetEntry(i);
        if i%2000 == 0: print "i = ",i
        
        for j in range(len(vars)):
            h_gg_s[j].Fill( getattr( tree_gg_s, vars[j] ) );

    for i in range( tree_pp_zp.GetEntries() ):
        tree_pp_zp.GetEntry(i);
        if i%2000 == 0: print "i = ",i

        for j in range(len(vars)):
            h_pp_zp[j].Fill( getattr( tree_pp_zp, vars[j] ) );
                
    ##### plot
    # scale and color
    for j in range(len(vars)):
        h_gg_s[j].Scale( 1./h_gg_s[j].Integral() );
        h_pp_zp[j].Scale( 1./h_pp_zp[j].Integral() );
        h_pp_zp[j].SetLineColor( 2 )
        h_pp_zp[j].SetLineWidth( 2 )
        h_gg_s[j].SetLineWidth( 2 )
        
    # on canvas
    for j in range(len(vars)):
        
        tmpcan = ROOT.TCanvas("can"+str(j),"can"+str(j),800,800);
        h_gg_s[j].SetMaximum( 1.2*max( h_gg_s[j].GetMaximum(), h_pp_zp[j].GetMaximum() ) );
        h_gg_s[j].SetMinimum( 0 );        
        h_gg_s[j].Draw();
        h_pp_zp[j].Draw("same");

        leg = ROOT.TLegend(0.25, 0.2, 0.5, 0.45)
        leg.SetBorderSize(0)
        leg.AddEntry(h_gg_s[j], "gg_s", "l")
        leg.AddEntry(h_pp_zp[j], "pp_zp", "l")
        leg.Draw()

        tmpcan.SaveAs("plots/"+vars[j]+".eps");
        tmpcan.SaveAs("plots/"+vars[j]+".png");        
        
        
        
        
        
        
        
        
