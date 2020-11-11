import ROOT

from random import randint,gauss
import time

def HistoProducerGaussOverExpo(HistoN,xbins,xmin,xmax,ybins,ymin,ymax,NB=10000,NS=1000,
                               expsets=[100.0,0.0,-0.0005,0.0002,100.0,0.0,-0.001,0.0005,100.0],
                               gausssets=[1.0,650,10,10.0,1.0,150,10,10.0]):
    Signal2DHistos=[]
    Backgr2DHistos=[]
    SigpBk2DHistos=[]
    Signal1DHistos=[]
    Backgr1DHistos=[]
    SigpBk1DHistos=[]
    Signal0DHistos=[]
    Backgr0DHistos=[]
    SigpBk0DHistos=[]
    for i in range(HistoN):
        BackgroundModel = ROOT.TH2F( 'bkg', 'My background model', xbins, xmin, xmax, ybins, ymin, ymax)
        F2D = ROOT.TF2("f2","([0]*expo(x,[1..2]))*([3]*expo(y,[4..5])+[6])",xmin, xmax, ymin, ymax)
        F2D.SetParameters(expsets[0],expsets[1],gauss(expsets[2],expsets[3]),expsets[4],expsets[5],gauss(expsets[6],expsets[7]),expsets[8])
        BackgroundModel.FillRandom("f2",NB)
        SignalModelforN = ROOT.TH2F( 'sigforN', 'My signal model', xbins, xmin, xmax, ybins, ymin, ymax)
        signalfforN = ROOT.TF2("signalf", "gausn(x, [0..2])*gausn(y, [3..5])",xmin, xmax, ymin, ymax)
        mux=gauss(gausssets[1],gausssets[2])
        muy=gauss(gausssets[5],gausssets[6])
        signalfforN.SetParameters(gausssets[0],mux,gausssets[3],gausssets[4],muy,gausssets[7])
        SignalModelforN.FillRandom("signalf",NS)
        BkgpSigModel = BackgroundModel.Clone("bkgpsig")
        BkgpSigModel.Sumw2()
        BkgpSigModel.Add(SignalModelforN,1.0)
        n2=[]; b2=[]; s2=[]
        n1=[]; b1=[]; s1=[]
        for i in range(1,BackgroundModel.GetNbinsX()):
            ni=0.0; si=0.0; bi=0.0
            for j in range(1,BackgroundModel.GetNbinsY()):
                n2.append(float(BkgpSigModel.GetBinContent(BkgpSigModel.GetBin(i,j))))
                s2.append(float(SignalModelforN.GetBinContent(SignalModelforN.GetBin(i,j))))
                b2.append(float(BackgroundModel.GetBinContent(BackgroundModel.GetBin(i,j))))
                ni+=float(BkgpSigModel.GetBinContent(BkgpSigModel.GetBin(i,j)))
                si+=float(SignalModelforN.GetBinContent(SignalModelforN.GetBin(i,j)))
                bi+=float(BackgroundModel.GetBinContent(BackgroundModel.GetBin(i,j)))
            n1.append(ni)
            s1.append(si)
            b1.append(bi)
        n0=[sum(n1)]
        s0=[sum(s1)]
        b0=[sum(b1)]
        Signal2DHistos.append(s2)
        Backgr2DHistos.append(b2)
        SigpBk2DHistos.append(n2)
        Signal1DHistos.append(s1)
        Backgr1DHistos.append(b1)
        SigpBk1DHistos.append(n1)
        Signal0DHistos.append(s0)
        Backgr0DHistos.append(b0)
        SigpBk0DHistos.append(n0)
    return Signal2DHistos, Backgr2DHistos, SigpBk2DHistos, Signal1DHistos, Backgr1DHistos, SigpBk1DHistos, Signal0DHistos, Backgr0DHistos, SigpBk0DHistos


def SimpleHistoProducerGaussOverExpo(HistoN,xbins,xmin,xmax,NB=10000,NS=1000,
                               expsets=[100.0,0.0,-0.0005],
                               gausssets=[1.0,650,10.0]):
    Signal1DHistos=[]
    Backgr1DHistos=[]
    SigpBk1DHistos=[]
    Signal0DHistos=[]
    Backgr0DHistos=[]
    SigpBk0DHistos=[]
    for i in range(HistoN):
        BackgroundModel = ROOT.TH1F( 'bkg', 'My background model', xbins, xmin, xmax)
        F1D = ROOT.TF1("f2","([0]*expo(x,[1..2]))",xmin, xmax)
        F1D.SetParameters(expsets[0],expsets[1],expsets[2])
        BackgroundModel.FillRandom("f2",NB)
        SignalModelforN = ROOT.TH1F( 'sigforN', 'My signal model', xbins, xmin, xmax)
        signalfforN = ROOT.TF1("signalf", "gausn(x, [0..2])",xmin, xmax)
        signalfforN.SetParameters(gausssets[0],gausssets[1],gausssets[2])
        SignalModelforN.FillRandom("signalf",NS)
        BkgpSigModel = BackgroundModel.Clone("bkgpsig")
        BkgpSigModel.Sumw2()
        BkgpSigModel.Add(SignalModelforN,1.0)
        n1=[]; b1=[]; s1=[]
        for i in range(1,BackgroundModel.GetNbinsX()):
            n1.append(float(BkgpSigModel.GetBinContent(BkgpSigModel.GetBin(i))))
            s1.append(float(SignalModelforN.GetBinContent(SignalModelforN.GetBin(i))))
            b1.append(float(BackgroundModel.GetBinContent(BackgroundModel.GetBin(i))))
        n0=[sum(n1)]
        s0=[sum(s1)]
        b0=[sum(b1)]
        Signal1DHistos.append(s1)
        Backgr1DHistos.append(b1)
        SigpBk1DHistos.append(n1)
        Signal0DHistos.append(s0)
        Backgr0DHistos.append(b0)
        SigpBk0DHistos.append(n0)
    return Signal1DHistos, Backgr1DHistos, SigpBk1DHistos, Signal0DHistos, Backgr0DHistos, SigpBk0DHistos
