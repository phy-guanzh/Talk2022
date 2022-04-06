import ROOT
import style
import CMS_lumi

ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
style.GoodStyle().cd()

colorful=[ROOT.kRed,ROOT.kBlue,ROOT.kGreen,ROOT.kCyan]
variables={
    'costheta1':{'name':'cos(#theta_{W_{1}})','rebin':1},
    'costheta2':{'name':'cos(#theta_{W_{2}})','rebin':1},
    'costheta11':{'name':'cos(#theta_{W_{1}})','rebin':1},
    'costheta22':{'name':'cos(#theta_{W_{2}})','rebin':1},
    'costheta13':{'name':'cos(#theta_{l1,W2})','rebin':1},
    'costheta23':{'name':'cos(#theta_{l2,W1})','rebin':1},
    'm_dphill_2':{'name':'#Delta#phi_{ll}','rebin':1},
    'm_met':{'name':'p_{T}^{miss} [GeV]','rebin':1},
    'm_met_2':{'name':'p_{T}^{miss} [GeV]','rebin':1},
    'm_mll':{'name':'m_{ll} (GeV)','rebin':1},
    'm_mll_2':{'name':'m_{ll} (GeV)','rebin':1},
    'm_mll_3':{'name':'m_{ll} (GeV)','rebin':1},
    'm_mlvlv':{'name':'m_{l#nul#nu}','rebin':1},
    'm_ptl1':{'name':'p_{T}^{l1} [GeV]','rebin':1},
    'm_ptl2':{'name':'p_{T}^{l2} [GeV]','rebin':1},
    'm_ptl2_2':{'name':'p_{T}^{l2} [GeV]','rebin':1},
    'etal1':{'name':'#eta_{l1}','rebin':1},
    'etal2':{'name':'#eta_{l2}','rebin':1},
    'm_ptg':{'name':'p_{T}^{photon} [GeV]','rebin':1},
}
frame=''
subdir='aqgc'
params = {
    'aqgc265':{'name':'polar','samples':['wwa_offshell_part'+frame,'wwa_offshell'+frame],'scale':[0.9722962041070317,0.9549387370405278]},
}

xoffsetstart = 0.0
yoffsetstart = 0.0
xoffset = 0.35
yoffset = 0.06
xpositions = [0.50,0.50,0.50,0.47,0.47,0.47,0.47,0.21,0.21,0.21,0.21]
ypositions = [0,1,2,0,1,2,3,0,1,2,3]
def draw_legend(x1,y1,hist,label,options):

    legend = ROOT.TLegend(x1+xoffsetstart,y1+yoffsetstart,x1+xoffsetstart + xoffset,y1+yoffsetstart + yoffset)

    legend.SetBorderSize(     0)
    legend.SetFillColor (     0)
    legend.SetTextAlign (    12)
    legend.SetTextFont  (    42)
    legend.SetTextSize  ( 0.040)
    print(label)
    legend.AddEntry(hist,label,options)

    legend.Draw("same")

    #otherwise the legend goes out of scope and is deleted once the function finishes
    hist.label = legend

def createPlot(ivar,i,name,samples,scale):
    c1 = ROOT.TCanvas("c1", "c1",5,50,500,500)
    #c1.SetLogy()
    f2=ROOT.TFile.Open('wwa_aqgc_ckm.root')
    print (ivar)
    h2=f2.Get(ivar+'_wwa_aqgc_ckm'+frame)
    print(ivar+'_wwa_aqgc_ckm'+frame,h2.Integral())
    h2.Sumw2()
    h2.SetLineColor(colorful[1])
    h2.SetMarkerColor(colorful[1])
    h2.SetLineWidth(3)


    f3=ROOT.TFile.Open('wwa_aqgc_ckm_centre.root')
    h3=f3.Get(ivar+'_wwa_aqgc_ckm_centre'+frame)
    h3.Sumw2()
    print(ivar+'_wwa_aqgc_ckm'+frame,h3.Integral())
    h3.SetLineColor(colorful[2])
    h3.SetMarkerColor(colorful[2])
    h3.SetLineWidth(3)
    

    h2.SetMaximum(max(h2.GetMaximum(),h3.GetMaximum())*1.4)

    plot_name={'wwa_offshell_part'+frame:'#bf{FT0=0}','wwa_offshell'+frame:'#bf{FT0=e-12}'}
    legend_count=0

    h2.SetMinimum(0)
    h2.SetTitle("")
    h2.GetYaxis().SetTitle("A.U.")
    h2.GetYaxis().SetTitleSize(0.04)
    h2.GetYaxis().SetTitleFont(62)
    h2.GetXaxis().SetTitleSize(0.05)
    h2.GetXaxis().SetTitleFont(62)
    h2.GetXaxis().SetTitle(variables[ivar]['name'])
    h2.Draw()
    draw_legend(xpositions[legend_count],0.84 - ypositions[legend_count]*yoffset,h2,plot_name['wwa_offshell_part'+frame],"l")
    legend_count+=1

    draw_legend(xpositions[legend_count],0.84 - ypositions[legend_count]*yoffset,h3,plot_name['wwa_offshell'+frame],"l")
    h3.Draw("same e")
    

    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "    Simulation"
    CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

    iPos = 0
    if( iPos==0 ): CMS_lumi.relPosX = 0.12
    iPeriod = 0
    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
    c1.RedrawAxis()
    c1.SaveAs(subdir+'/'+ivar+'_'+i+frame+'_shape.png')
    c1.SaveAs(subdir+'/'+ivar+'_'+i+frame+'_shape.pdf')
    c1.Close()
    f2.Close()
    
for ivar in variables:
    for i in params:
        createPlot(ivar,i,params[i]['name'],params[i]['samples'],params[i]['scale'])

'''
for isample in range(0,len(samples)):
    for ivar in variables:
        c1 = ROOT.TCanvas("c1", "c1",5,50,500,500)
        c1.SetLogy()
        tl = ROOT.TLegend(0.65,0.78,0.87,0.87)
        h1=f1.Get(ivar+'_'+samples[isample])
        h2=f2.Get(ivar+'_'+samples2[isample])
        #h1=f1.Get('m_ptl1_'+samples[isample])
        #h2=f2.Get('m_ptl2_'+samples2[isample])
        norm = h2.Integral()
        #h1.Scale(h2.Integral()/h1.Integral())
        #h1.Rebin(2)
        #h2.Rebin(2)
        if h1.GetMaximum() < h2.GetMaximum():
            h1.SetMaximum(h2.GetMaximum()*1.2)
        else:
            h1.SetMaximum(h1.GetMaximum()*1.2)
        h1.SetLineColor(ROOT.kRed)
        h1.SetMarkerColor(ROOT.kRed)
        h1.GetXaxis().SetTitle('p_{T} [GeV]')#variables[ivar]
        h1.SetTitle('VBS W^{-}W^{-}: pt')#+variables[ivar]
        h2.SetLineColor(ROOT.kBlue)
        h2.SetMarkerColor(ROOT.kBlue)
        h1.Draw("")
        h2.Draw("same E")

        tl.AddEntry(h1,'electron','l')
        tl.AddEntry(h2,'muon','l')
        tl.Draw('same')
        c1.SaveAs(ivar+'_'+samples[isample]+".pdf")
        c1.Close()
'''
