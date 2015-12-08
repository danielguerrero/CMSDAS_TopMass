#!/usr/bin/env python
import ROOT
from ROOT import *
import os,sys,pprint

def doPeakFit(h=None,minToFit=None,maxToFit=None,outputName=None):

    # Set the stats off 
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)
    gStyle.SetTickLength(0.03)

    #Get the log(E) histogram and change range around the peak
    hFit = h.Clone()
    hFit.GetXaxis().SetRangeUser(minToFit,maxToFit)     

    #Define the fit function and make the fit to log(E)
    gauss = TF1("gauss", "gaus", minToFit, maxToFit)
    gauss.SetLineColor(kBlue)
    gauss.SetLineWidth(3)
    gauss.SetLineStyle(1)
    hFit.Fit("gauss","EQ") 

    #Make a Pull histogram for your fit    
    hPull = h.Clone()
    hPull.GetXaxis().SetRangeUser(minToFit,maxToFit)
    for bin in range(0, hFit.GetNbinsX()):
       hPull.SetBinContent(bin, 0.)
       hPull.SetBinError(bin, 0.)
    for ibin in range(0, hFit.GetNbinsX()):
       binCont = hFit.GetBinContent(ibin)
       binErr = hFit.GetBinError(ibin)
       valIntegral = gauss.Eval( hFit.GetBinCenter(ibin))
       if binErr !=0:
         pull = (binCont-valIntegral)/binErr
         hPull.SetBinContent(ibin, pull)
         hPull.SetBinError(ibin, 1)

    #Create a canvas for plotting your histograms
    c=TCanvas('c','c')
    p1=ROOT.TPad('p1','p1',0.,0.3,1.0,1.0)
    p2=ROOT.TPad('p2','p2',0.,0.,1.0,0.3)
    p1.Draw()    
    p2.Draw()

    #Edit the pad of the fit
    p1.cd()
    p1.SetBorderMode(0)
    p1.SetBorderSize(2)
    p1.SetTickx(1)
    p1.SetTicky(1)
    p1.SetTopMargin(0.13)
    p1.SetBottomMargin(0.02)
    gauss.SetLineWidth(3)
    gauss.SetLineColor(kBlue)
    gauss.SetLineStyle(1)
    hFit.SetMarkerStyle(8)
    hFit.GetYaxis().SetTitleSize(0.062)
    hFit.GetYaxis().SetLabelSize(0.062)
    hFit.GetYaxis().SetTitleOffset(0.62)
    hFit.GetXaxis().SetLabelOffset(1)
    hFit.GetYaxis().SetTitle("1/E dN_{bjets}/dlog(E)")
    hFit.GetXaxis().SetTitle("log(E)")
    hFit.SetLineColor(kBlack)
    hPull.SetLineColor(kBlack)
    hFit.SetMarkerColor(kBlack)
    hPull.SetMarkerColor(kBlack)
    hFit.Draw()

    #Get Fit Parameters
    mean=gauss.GetParameter(1)
    meanErr=gauss.GetParError(1)
    sigma=gauss.GetParameter(2)
    sigmaErr=gauss.GetParError(2)
    chi2=gauss.GetChisquare()
    NDF=gauss.GetNDF()
    chi2ndf=chi2/NDF
    #Calculate the uncalibrated Energy peak position and its uncertainty
    Ereco=exp(mean)
    Err=abs(Ereco*meanErr)
    Eup=Ereco+Err
    #Calculate the uncalibrated top mass and its uncertainty
    mtop=Ereco+sqrt((80.4*80.4)-(4.8*4.8)+(Ereco*Ereco))
    mtopup=Eup+sqrt((80.4*80.4)-(4.8*4.8)+(Eup*Eup))
    mtopErr=abs(mtop-mtopup)

    #Create some labels about the statistics
    caption1=TLatex()
    caption1.SetTextSize(0.038)
    caption1.SetTextFont(42)
    caption1.SetNDC()
    caption1.DrawLatex(0.75,0.82,'Fit Results')
    caption1.DrawLatex(0.7,0.78,'Mean=%3.4f #pm %3.4f'%(mean,meanErr))
    caption1.DrawLatex(0.7,0.74,'Width=%3.4f #pm %3.4f '%(sigma,sigmaErr))
    caption1.DrawLatex(0.73,0.69,'#chi^{2}/ndf=%3.4f'%(chi2ndf))
    caption2=TLatex()
    caption2.SetTextSize(0.05)
    caption2.SetTextFont(42)
    caption2.SetNDC()  
    caption2.DrawLatex(0.40,0.44,'Uncalibrated Measurement')
    caption2.DrawLatex(0.38,0.39,'E_{peak} = %3.4f #pm %3.4f GeV'%(Ereco,Err))
    caption2.DrawLatex(0.39,0.33,'m_{t} = %3.4f #pm %3.4f GeV'%(mtop,mtopErr))

    #CMS labels
    tlat1 = TLatex()
    tlat1.SetNDC()
    tlat1.SetTextFont(60)
    tlat1.SetTextSize(0.10)
    tlat1.SetTextAlign(31)
    prelim_text1 = 'CMS'
    tlat1.DrawLatex(0.25, 0.77, prelim_text1)
    tlat1 = TLatex()
    tlat1.SetNDC()
    tlat1.SetTextFont(42)
    tlat1.SetTextSize(0.085)
    tlat1.SetTextAlign(31)
    prelim_text2 ='#it{Simulation}'
    tlat1.DrawLatex(0.30, 0.70, prelim_text2)
    prelim_text3 =' 8 TeV'
    tlat1.DrawLatex(0.90, 0.90, prelim_text3)

    #Edit the pad for the pull
    p2.cd()
    p2.SetGridy()
    p2.SetBorderMode(0)
    p2.SetBorderSize(2)
    p2.SetTickx(1)
    p2.SetTicky(1)
    p2.SetTopMargin(0.05)
    p2.SetBottomMargin(0.3)
    hPull.SetMarkerStyle(8)
    hPull.GetYaxis().SetTitle("#frac{Data-Fit}{Uncertainty}")
    hPull.GetYaxis().SetTitleSize(0.140)
    hPull.GetYaxis().SetLabelSize(0.140)
    hPull.GetXaxis().SetTitleSize(0.160)
    hPull.GetXaxis().SetLabelSize(0.150)
    hPull.GetYaxis().SetTitleOffset(0.23)
    hPull.GetXaxis().SetTitleOffset(0.8)
    hPull.Draw("e")
    hPull.GetXaxis().SetTitle("log(E)")

    #save and delete
    c.SaveAs( outputName + '.pdf' )
    c.SaveAs( outputName + '.C' )    
    del c
    gauss.IsA().Destructor(gauss)
    del caption1,caption2

    #all done here ;)
    return mean,meanErr,sigma,sigmaErr,Ereco,Err

def getPeakValue(h=None):
           histoname="Eb_emu_rec_fit"
           #Get the bin of the mean in the energy spectrum
           m=h.GetMean()
           rms=h.GetRMS()
           i=h.GetXaxis().FindBin(m)
           imean=h.GetBinCenter(i)
           #Get the range where you are going to make the initial fit
           isigma=rms
           minlimit = imean - isigma
           maxlimit = imean + isigma
           #Apply an initial fit in order to locate the peak as a first approach
           mean1,meanUnc1,sigma1,sigma1Unc,Ereco1,Err1=doPeakFit(h=h,minToFit=minlimit,maxToFit=maxlimit,outputName=histoname)
           j=h.GetXaxis().FindBin(mean1)
           jmean=h.GetBinCenter(j)
           newmean = jmean     
           newsigma = isigma
           #Redefine your range of location of the peak for the final fit
           minlimit = newmean - newsigma
           maxlimit = newmean + newsigma
           #Make the refining fit  
           mean,meanUnc,sigma,sigmaUnc,Ereco,Err=doPeakFit(h=h,minToFit=minlimit,maxToFit=maxlimit,outputName=histoname)      
           return Ereco,Err

def Histo_PEs(histopes=None,rootname=None,xlabel=None,ylabel=None):
    canvas=TCanvas("canvas","canvas")
    gStyle.SetOptStat(0)  
    histopes.SetMarkerStyle(8)
    histopes.GetYaxis().SetTitleSize(0.047)
    histopes.GetYaxis().SetLabelSize(0.047)
    histopes.GetXaxis().SetTitle("%s"%(xlabel))           
    histopes.GetYaxis().SetTitle("%s"%(ylabel))
    histopes.Draw()
    histopes.SaveAs("%s.root"%(rootname))
    canvas.SaveAs("%s.pdf"%(rootname))
    del canvas

def main(argv=None):

           # Get the variables = rootfile name and histogram name
           filename1 = sys.argv[1]
           h1 = sys.argv[2]   
           print "... processing", filename1
           from os import path
           if not path.isfile(filename1):
               print "Help, file doesnt exist"
               exit(-1)

           #Get the histograms in the template sample 
           f = ROOT.TFile(filename1, "read")
           eblogHisto = f.Get(h1)
           #Book histograms
           hist = TH1D("hist","hist",80,3.4,7.4)
           Emeasured = TH1D("Emeasured","Emeasured",1000,40.,90.)
           Uncertainty = TH1D("Uncertainty","Uncertainty", 400, 0.,2.)
           hist.Sumw2()
           #Calculate the energy peak position in the big MC sample
           Ebreco,Ebrr=getPeakValue(h=eblogHisto)
           print Ebreco,Ebrr

           # Create the 10000 pseudoexperiments, make the measurement and fill distributions
           r=TRandom3()
           r.SetSeed(1)
           #Run a loop for generating the pseudoexperiments
           for num in range (0,10000):
              hist.Reset('ICE') 
              #Produce the energy spectrum in hte PE
              nbin=eblogHisto.GetNbinsX()             
              for j in range(0,nbin+1):
                   E = exp(eblogHisto.GetBinCenter(j))
                   fluct = r.Poisson(E*eblogHisto.GetBinContent(j))
                   hist.SetBinError(j, sqrt(fluct)/E)                
                   hist.SetBinContent(j, fluct /E )                                 
              #Get the information about the fit in the PE
              Ereco,Err=getPeakValue(h=hist)
              #Fill the histograms that you are interested          
              Emeasured.Fill(Ereco)
              Uncertainty.Fill(Err)  
              print "Number of pseudoexperiment:", num
           # Create a canvas and plot the histograms of the pseudo-experiments result           
           Histo_PEs(histopes=Emeasured,rootname="Emeasured",xlabel="Energy [GeV|]", ylabel="Pseudoexperiments")
           Histo_PEs(histopes=Uncertainty,rootname="Uncertainty",xlabel="Uncertainty [GeV]", ylabel="Pseudoexperiments")
               
if __name__ == "__main__":
    sys.exit(main())
