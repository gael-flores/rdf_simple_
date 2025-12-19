import ROOT
import sys
from array import array
import numpy as np
import pickle
import ctypes
import math
import common.tdrstyle as tdrstyle
import common.CMS_lumi as CMS_lumi
from python.VHTools.py_helpers import *
sty = tdrstyle.setTDRStyle()
from scipy.stats import chi2
### = commented out to test tdrStyle, uncomment if not using
import matplotlib.pyplot as plt
import mplhep as mh

drawprelim = True

@contextmanager
def open_root_file(filename, mode="READ"):
    """
    Context manager for safely opening and closing a ROOT TFile.

    This ensures that the ROOT file is automatically closed when the
    context exits, even in the event of an exception.

    Parameters:
    ----------
    filename : str
        Path to the ROOT file to open.

    mode : str, optional
        Mode in which to open the file (default is "READ").
        Typical options include:
        - "READ"   : Open existing file for reading
        - "UPDATE" : Open existing file for reading/writing
        - "RECREATE": Create new file, overwriting if exists
        - "NEW"    : Create new file, fail if exists

    Yields:
    ------
    ROOT.TFile
        The opened ROOT file object.

    Example:
    --------
    >>> with open_root_file("data.root") as f:
    ...     h = f.Get("myHistogram")
    ...     h.Draw()
    """
    file = ROOT.TFile.Open(filename, mode)
    try:
        yield file
    finally:
        if file and file.IsOpen():
            file.Close()

class plotter_base(object):

    def __init__(self):
        self.fillstyle=1001
        self.linestyle=1
        self.linecolor=1
        self.linewidth=2
        self.fillcolor=ROOT.kOrange-3
        self.markerstyle=20
        self.corrFactors=[]

    def addCorrectionFactor(self,value,model):
        corr=dict()
        corr['value']=value
        corr['model']=model
        self.corrFactors.append(corr)

    def setLineProperties(self,linestyle,linecolor,linewidth):
        self.linestyle=linestyle
        self.linecolor=linecolor
        self.linewidth=linewidth 

    def setFillProperties(self,fillstyle,fillcolor):
        self.fillstyle=fillstyle
        self.fillcolor=fillcolor

    def setMarkerProperties(self,markerstyle):
        self.markerstyle=markerstyle
        
        
    def array1d(self,var,cuts,model,include_overflow=False):
        hist=self.hist1d(var,cuts,model,titlex = "",units = "")
        axis=hist.GetXaxis()
        num_bins = hist.GetNbinsX()+2
        data = np.ndarray(num_bins, dtype=np.float64, buffer=hist.GetArray())
        w2   = np.ndarray(num_bins, dtype=np.float64, buffer=hist.GetSumw2().GetArray())
        edges = np.array([axis.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX()+1)])
        edges=np.append(edges,axis.GetBinUpEdge(num_bins))
        if include_overflow==False:
            return edges.copy(),data[1:-1].copy(),w2[1:-1].copy()
        else:
            return edges.copy(),data.copy(),w2.copy()
    def array2d(self,var1,var2,cuts,model,include_overflow=False):
        h=self.hist2d(var1,var2,cuts,model,titlex = "",unitsx = "",titley="",unitsy="")        
        nx = h.GetNbinsX() + 2
        ny = h.GetNbinsY() + 2
        data = np.ndarray((nx, ny), buffer=h.GetArray(), dtype=np.float64)
        w2   = np.ndarray((nx, ny), buffer=h.GetSumw2().GetArray(), dtype=np.float64)
        xaxis=h.GetXaxis()
        yaxis=h.GetYaxis()
        xedges = np.array([xaxis.GetBinLowEdge(i) for i in range(1, h.GetNbinsX() + 1)])
        xedges=np.append(xedges,xaxis.GetBinUpEdge(nx))                     
        yedges = np.array([yaxis.GetBinLowEdge(i) for i in range(1, h.GetNbinsY() + 1)])
        yedges=np.append(yedges,yaxis.GetBinUpEdge(ny))
        if include_overflow==False:
            return xedges.copy(),yedges.copy(),data[1:-1.1:-1].copy(),w2[1:-1,1:-1].copy()
        else:
            return xedges.copy(),yedges.copy(),data.copy(),w2.copy()        




class rdf_plotter(plotter_base):
    def __init__(self,file,isMC=False,weight = "1.0",tree='Events',defaultCuts = "1.0",report='report'):
        self.rdf = ROOT.RDataFrame(tree,file)
        self.weight=weight
        self.isMC=isMC
        self.defaultCuts = defaultCuts
        self.file = file
        self.report = report
        super(rdf_plotter,self).__init__()
        #if MC read the sum of weights and weigh the event
        if isMC:
            with open_root_file(file) as f:
                t = f.Get('Runs')
                sumw = 0.0
                for event in t:
                    sumw += event.genEventSumw
                self.weight = self.weight+'*(genWeight/{})*(sample_sigma)*(Pileup_weight)'.format(sumw if sumw>0 else 1.0)

    def readReport(self) -> dict[str]:
        out = {}
        with open_root_file(self.file) as f:
            report = f.Get(self.report)
            if not report:
                raise RuntimeError(f"Report {self.report} not found in {self.file}")
            
            cuts = [b.GetName() for b in report.GetListOfBranches()]
            out = {c: 0 for c in cuts}
            for event in report:
                for c in cuts:
                    out[c] += getattr(event, c, 0)

            return out
        

    def define(self, var, definition):
        if var in self.rdf.GetColumnNames():
            self.rdf = self.rdf.Redefine(var, definition)           
        else:    
            self.rdf = self.rdf.Define(var, definition)

    def redefine(self, var, definition):
        self.rdf = self.rdf.Redefine(var, definition)
    
    def display(self,var):
        self.rdf.Display(var).Print()

    def filter(self,condition):
        '''
        Inputs
        ======
        condition(string): condition to be used for filtering the dataframe that evaluates to True (1) or False (0)
        new_var(string): name of the new column  
        
        Returns
        =======
        RDataFrame with filter applied
        '''
        self.rdf = self.rdf.Filter(condition)

        
    def hist1d(self,var,cuts,model,titlex = "",units = ""):
        corrString="1.0"
        for corr in self.corrFactors:
            corrString = corrString+"*("+str(corr['value'])+")" 
        c = "("+self.defaultCuts+")*("+cuts+")*"+self.weight+"*("+corrString+")"
        rdf=self.rdf.Define('plot_weight',c)
        h=rdf.Histo1D(model,var,'plot_weight')
        h.SetLineStyle(self.linestyle)
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        if units=="":
            h.GetXaxis().SetTitle(titlex)
        else:
            h.GetXaxis().SetTitle(titlex+ " ["+units+"]")
        for n in range(1, h.GetNbinsX()+1):
            if h.GetBinContent(n) < 0:
                h.SetBinContent(n, 0)
        return h.GetPtr()



    
    def profile1d(self,var1,var2,cuts,model,titlex = "",unitsx = "",titley = "",unitsy = ""):
        corrString="1.0"
        for corr in self.corrFactors:
            corrString = corrString+"*("+str(corr['value'])+")" 
        c = "("+self.defaultCuts+")*("+cuts+")*"+self.weight+"*("+corrString+")"
        rdf=self.rdf.Define('plot_weight',c)
        h=rdf.Profile1D(model,var1,var2,'plot_weight')
        h.SetLineStyle(self.linestyle)
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        if unitsx=="":
            h.GetXaxis().SetTitle(titlex)
        else:
            h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        if unitsy=="":
            h.GetYaxis().SetTitle(titley)
        else:
            h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        return h.GetPtr()


    def hist2d(self,var1,var2,cuts,model,titlex = "",unitsx = "",titley="",unitsy=""):
        corrString="1.0"
        for corr in self.corrFactors:
            corrString = corrString+"*("+str(corr['value'])+")" 
        c = "("+self.defaultCuts+")*("+cuts+")*"+self.weight+"*("+corrString+")"
        rdf=self.rdf.Define('plot_weight',c)
        h=rdf.Histo2D(model,var1,var2,'plot_weight')

        #h.Sumw2(0)
        #h.SetBinErrorOption(ROOT.TH1D.kPoisson)
        h.SetLineStyle(self.linestyle)
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        if unitsx=="":
            h.GetXaxis().SetTitle(titlex)
        else:
            h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        if unitsy=="":
            h.GetYaxis().SetTitle(titley)
        else:
            h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        return h.GetPtr()

    def unrolled2d(self,var1,var2,cuts,model):
        corrString="1.0"
        for corr in self.corrFactors:
            corrString = corrString+"*("+str(corr['value'])+")" 
        c = "("+self.defaultCuts+")*("+cuts+")*"+self.weight+"*("+corrString+")"
        rdf=self.rdf.Define('plot_weight',c)
        h=rdf.Histo2D(model,var1,var2,'plot_weight')

        hu=ROOT.TH1D(h.GetName(),h.GetTitle(),h.GetNbinsX()*h.GetNbinsY(),0,h.GetNbinsX()*h.GetNbinsY())
        for i in range(1,h.GetNbinsX()+1):
            for j in range(1,h.GetNbinsY()+1):
                b=h.GetBin(i,j)
                bin1d=(i-1)*h.GetNbinsX()+j
                hu.SetBinContent(bin1d,h.GetBinContent(b))
                hu.SetBinError(bin1d,h.GetBinError(b))
            
        hu.SetLineStyle(self.linestyle)
        hu.SetLineColor(self.linecolor)
        hu.SetLineWidth(self.linewidth)
        hu.SetFillStyle(self.fillstyle)
        hu.SetFillColor(self.fillcolor)
        hu.SetMarkerStyle(self.markerstyle)
        return hu


    

class merged_plotter(plotter_base):

    def __init__(self, plotters):
        super(merged_plotter,self).__init__()
        self.plotters = plotters

    def define(self, var, definition):
        for plotter in self.plotters:
            plotter.define(var, definition)

    def filter(self, condition):
        for plotter in self.plotters:
            plotter.filter(condition)

    def redefine(self, var, definition):
        for plotter in self.plotters:
            plotter.redefine(var, definition)

    def hist1d(self,var,cuts,model,titlex = "",units = ""):
        h = None
        for plotter in self.plotters:
            if h is None:
                h = plotter.hist1d(var, cuts, model, titlex, units)
            else:
                h.Add(plotter.hist1d(var, cuts, model, titlex, units))
        if h is None:
            tmprdf = ROOT.RDataFrame(1)
            tmprdf = tmprdf.Define("x", "0")
            h = tmprdf.Define("goodX", "x!=0").Histo1D(model, "goodX")
        h.SetLineStyle(self.linestyle)
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        if units=="":
            h.GetXaxis().SetTitle(titlex)
        else:
            h.GetXaxis().SetTitle(titlex+ " ["+units+"]")
        return h


    def hist2d(self,var1,var2,cuts,model,titlex = "",unitsx = "",titley="",unitsy=""):
        h = None
        for plotter in self.plotters:
            if h is None:
                h = plotter.hist2d(var1,var2, cuts, model, titlex, unitsx, titley, unitsy)
            else:
                h.Add(plotter.hist2d(var1,var2, cuts, model, titlex, unitsx, titley, unitsy))
        if h is None:
            return h
        else:
            h.SetLineStyle(self.linestyle)
            h.SetLineColor(self.linecolor)
            h.SetLineWidth(self.linewidth)
            h.SetFillStyle(self.fillstyle)
            h.SetFillColor(self.fillcolor)
            h.SetMarkerStyle(self.markerstyle)
            if unitsx=="":
                h.GetXaxis().SetTitle(titlex)
            else:
                h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
            if unitsy=="":
                h.GetYaxis().SetTitle(titley)
            else:
                h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
            return h

    def unrolled2d(self,var1,var2,cuts,model):
        h = None
        for plotter in self.plotters:
            if h is None:
                h = plotter.unrolled2d(var1,var2, cuts, model)
            else:
                h.Add(plotter.unrolled2d(var1,var2, cuts, model))
        if h is None:
            return h
        else:
            h.SetLineStyle(self.linestyle)
            h.SetLineColor(self.linecolor)
            h.SetLineWidth(self.linewidth)
            h.SetFillStyle(self.fillstyle)
            h.SetFillColor(self.fillcolor)
            h.SetMarkerStyle(self.markerstyle)
        return h
                               
    def getArray(self, var, cuts):
        """
        Applies cuts to RDF and extracts numpy array of values for a specific variable.

        Parameters:
            vars (str): variable
            cuts (str): Selection cuts

            range (tuple) : range over

        Returns:
            list of values corresponding to cuts
        """        
        # Initialize empty lists to store the data
        combined_var = []

        # Iterate over the plotters and collect data
        for plotter in self.plotters:
            # Define weight expression
            corrString = "1.0"
            for corr in plotter.corrFactors:
                corrString += f"*({corr['value']})"

            # Apply cuts and filters based on the provided cuts and lumi
            weight_expr = f"({plotter.defaultCuts})*({cuts})*{plotter.weight}*({corrString})"
            rdf = plotter.rdf.Define("plot_weight", weight_expr)
            rdf = rdf.Filter("plot_weight")

            # Get data from the current plotter
            numpy_var = rdf.AsNumpy([var])[var]

            # Collect the data
            combined_var.extend(numpy_var)

        # If no data points were collected, return None
        if len(combined_var) == 0:
            return None

        

        return combined_var
    
    def display(self,var):
        for plotter in self.plotters:
            plotter.display(var)

    def readReport(self):
        out = {}
        for plotter in self.plotters:
            r = plotter.readReport()
            for cut in r:
                if cut in out:
                    out[cut] += r[cut]
                else:
                    out[cut] = r[cut]
        return out



#Tricky, use with care!
class background_plotter(merged_plotter):
    def __init__(self,cutsSR,cutsCR,cutsSSB,cutsCSB,plotters):
        self.cutsSR=cutsSR
        self.cutsCR=cutsCR
        self.cutsSSB=cutsSSB
        self.cutsCSB=cutsCSB
        super(background_plotter,self).__init__(plotters)
        super(background_plotter,self).define('scaleVar','1.0')

    def getScale(self,cuts):       
        numeratorCuts = cuts.replace(self.cutsSR,self.cutsSSB)
        denominatorCuts = cuts.replace(self.cutsSR,self.cutsCSB)
        print(numeratorCuts)
        print(denominatorCuts)
        
        hNum=super(background_plotter,self).hist1d('scaleVar',numeratorCuts,('a','a',10,-5,5))
        #replace cuts with cuts from control region
        hDenom=super(background_plotter,self).hist1d('scaleVar',denominatorCuts,('a','a',10,-5,5))
        if hDenom.Integral()==0:
            print("Error , denominator has zero events, returning scale=1")
            return 1.0
        else:
            return hNum.Integral()/hDenom.Integral()        
        
    def hist1d(self,var,cuts,model,titlex = "",units = ""):
        h=super(background_plotter,self).hist1d(var,cuts.replace(self.cutsSR,self.cutsCR),model,titlex,units)
        h.Scale(self.getScale(cuts))
        return h
    
    def hist2d(self,var1,var2,cuts,model,titlex = "",unitsx = "",titley="",unitsy=""):
        h=super(background_plotter,self).hist2d(var1,var2,cuts.replace(self.cutsSR,self.cutsCR),model,titlex,unitsx,titley,unitsy)
        h.Scale(self.getScale(cuts))
        return h

    def unrolled2d(self,var1,var2,cuts,model):
        h=super(background_plotter,self).unrolled2d(var1,var2,cuts.replace(self.cutsSR,self.cutsCR),model)      
        h.Scale(self.getScale(cuts))
        return h
    
class stack_plotter(object):
    def __init__(self,mode='stack',defaultCut="1"):
        self.plotters = []        
        self.log=False
        self.defaultCut=defaultCut
        self.mode=mode
    def setLog(self,doLog):
        self.log=doLog

    def add_plotter(self,plotter,name="",label = "label",typeP = "background"):
        packet = {'plotter':plotter,
                  'type':typeP,
                  'label':label,
                  'name':name}
        self.plotters.append(packet)
        
    def define(self, var, definition):
        for plotter in self.plotters:
            plotter['plotter'].define(var, definition)

    def redefine(self, var, definition):
        for plotter in self.plotters:
            plotter['plotter'].redefine(var, definition)

            
    def __getattr__(self, name):

        def fill_histograms(*args,**kwargs):    
            bkg_histograms=[]
            sig_histograms=[]
            data_histograms=[]
            
        
            #extract the histograms and build the stack
            for p in self.plotters:
                if hasattr(p['plotter'], name):
                    method = getattr(p['plotter'], name)
                    if p['type']=='background':
                        bkg_histograms.append(method(*args,**kwargs))
                        bkg_histograms[-1].SetName(p['name'])
                    elif p['type']=='signal':
                        sig_histograms.append(method(*args,**kwargs))
                        sig_histograms[-1].SetName(p['name'])
                                                   
                    elif p['type']=='data':
                        data_histograms.append(method(*args,**kwargs))
                        data_histograms[-1].SetName(p['name'])

            maxYields = [s.GetMaximum() for s in sig_histograms+bkg_histograms+data_histograms]
            integrals = [s.GetMaximum()/s.Integral() for s in sig_histograms+bkg_histograms+data_histograms]
            maxYield = max(maxYields)*1.2
            maxNormYield = max(integrals)*1.5
            #ok now stack all background hists if stack is the mode
            if self.mode=='stack':
                for i in range(0,len(bkg_histograms)-1):
                    for j in range(i,len(bkg_histograms)):
                        bkg_histograms[i].Add(bkg_histograms[j])
            #also add all signals separately to the background stack
            if len(bkg_histograms)>0:
                for h in sig_histograms:
                    h.Add(bkg_histograms[0])
                
                
        
            canvas = ROOT.TCanvas("canvas", "", 800, 600)
            canvas.SetTopMargin(0.055) # testing tdr style
            canvas.SetBottomMargin(0.12) # testing
            canvas.SetLeftMargin(0.12) # testing

            canvas.cd()
            #first the background
            for i in range(0,len(bkg_histograms)):
                if i==0:
                    if self.mode=='stack':
                        bkg_histograms[i].GetYaxis().SetRangeUser(0,maxYield)
                        bkg_histograms[i].Draw("HIST")
                        
                    else:
                        bkg_histograms[i].GetYaxis().SetRangeUser(0,maxNormYield)                        
                        bkg_histograms[i].DrawNormalized("HIST")                    
                else:    
                    if self.mode=='stack':
                        bkg_histograms[i].Draw("HIST,SAME")
                    else:
                        bkg_histograms[i].DrawNormalized("HIST,SAME")                    
            #then the signals
            for i in range(0,len(sig_histograms)):
                if i==0 and len(bkg_histograms)==0:
                    if self.mode=='stack':
                        sig_histograms[i].GetYaxis().SetRangeUser(0,maxYield)
                        sig_histograms[i].Draw("HIST")
                    else:
                        sig_histograms[i].GetYaxis().SetRangeUser(0,maxNormYield)
                        sig_histograms[i].DrawNormalized("HIST")                        
                else:
                    if self.mode=='stack':
                        sig_histograms[i].Draw("HIST,SAME")
                    else:
                        sig_histograms[i].DrawNormalized("HIST,SAME")                    
            
            #then the data
            for i in range(0,len(data_histograms)):
                if i==0 and len(bkg_histograms)==0 and len(sig_histograms)==0:
                    if self.mode=='stack':
                        data_histograms[i].GetYaxis().SetRangeUser(0,maxYield)
                        data_histograms[i].Draw("")
                    else:
                        data_histograms[i].GetYaxis().SetRangeUser(0,maxNormYield)
                        data_histograms[i].DrawNormalized("")
                else:        
                    if self.mode=='stack':
                        data_histograms[i].Draw("SAME")
                    else:
                        data_histograms[i].DrawNormalized("SAME")

            #Then the legend            
            legend = ROOT.TLegend(0.61,0.64,0.91,0.94,"","brNDC") # if using tdrstyle
            legend.SetBorderSize(0)
            legend.SetLineColor(1)
            legend.SetLineStyle(1)
            legend.SetLineWidth(1)
            legend.SetFillColor(0)
            legend.SetFillStyle(0)
            legend.SetTextFont(42)

            legend.SetFillColor(ROOT.kWhite)
            i=0
            for p in self.plotters:
                if p['type']=='data':
                    legend.AddEntry(data_histograms[i], p['label'], 'p')
                    i=i+1
            i=0
            for p in self.plotters:
                if p['type']=='signal':
                    legend.AddEntry(sig_histograms[i], p['label'], 'l')
                    i=i+1

            i=0
            for p in self.plotters:
                if p['type']=='background':
                    legend.AddEntry(bkg_histograms[i], p['label'], 'lf')
                    i=i+1
                
            legend.SetFillStyle(0)
            legend.SetBorderSize(0)
            ROOT.SetOwnership(legend, False)
            legend.Draw()
            canvas.Update()
            
            return {'canvas': canvas, 'legend': legend,'hists': sig_histograms+bkg_histograms+data_histograms}
    
        return fill_histograms



#MPLHEP specific plotter
class mplhep_plotter(object):
    def __init__(self,label='Preliminary',lumi=137.62,defaultCut="1",stack=True):
        self.plotters = []        
        self.defaultCut=defaultCut
        mh.style.use('CMS')
        self.label=label
        self.lumi=lumi
        self.stack=stack
        
    def add_plotter(self,plotter,name='name',label = "label",typeP = "background",color='black'):
        packet = {'plotter':plotter,
                  'type':typeP,
                  'label':label,
                  'name':name,
                  'color':color}
        self.plotters.append(packet)
        
    def define(self, var, definition):
        for plotter in self.plotters:
            plotter['plotter'].define(var, definition)

    def redefine(self, var, definition):
        for plotter in self.plotters:
            plotter['plotter'].redefine(var, definition)

    def plot1d(self,var,cuts,model,alpha=1.0,xlabel="",xunits="",legend_loc='upper right',show=True):
        background_hists=[]
        background_edges=[]
        background_w2=[]
        background_labels=[]
        background_colors=[]
        data_hists=[]
        data_edges=[]
        data_w2=[]
        data_labels=[]
        data_colors=[]
        signal_hists=[]
        signal_edges=[]
        signal_w2=[]
        signal_labels=[]
        signal_colors=[]

        bkgExists=False        
        for p in self.plotters:
            if p['type']=='data':
                edges,data,w2=p['plotter'].array1d(var,cuts,model)
                data_hists.append(data)
                data_edges.append(edges)
                data_w2.append(w2)
                data_labels.append(p['label'])
                data_colors.append(p['color'])                                    
            elif p['type']=='background':
                edges,data,w2=p['plotter'].array1d(var,cuts,model)
                if bkgExists==False:
                    background_sum=data
                    background_sumw2=w2
                    bkgExists=True
                else:
                    background_sum=background_sum+data
                    background_sumw2=background_sumw2+w2

                background_hists.append(data)
                background_edges.append(edges)
                background_w2.append(w2)
                background_labels.append(p['label'])
                background_colors.append(p['color'])                                                    
        for p in self.plotters:                   
            if p['type']=='signal':
                edges,data,w2=p['plotter'].array1d(var,cuts,model)
                if self.stack==True:
                    signal_hists.append(data+background_sum)
                    signal_w2.append(w2+background_sumw2)
                else:
                    signal_hists.append(data)
                    signal_w2.append(w2)                    
                signal_edges.append(edges)
                signal_labels.append(p['label'])
                signal_colors.append(p['color'])                                    


                
        fig,ax = plt.subplots()
        if len(signal_hists)>0:
            mh.histplot(signal_hists,background_edges[0],
                        histtype='step',
                        stack=False,
                        label=signal_labels,
                        color=signal_colors,
                        yerr=None,
                        sort='label',
                        ax=ax
                        )                   
        if len(background_hists)>0:
            #plot background stack            
            mh.histplot(background_hists,background_edges[0],
                        histtype='fill',
                        alpha = 1.0 if self.stack==True else 0.7,
                        stack=self.stack,
                        label=background_labels,
                        sort='label',
                        ax=ax
                        )
            #plot background error band in a custom way (since we use old version of mplhep)
            if self.stack:
                ax.fill_between(background_edges[0],
                                np.append(background_sum-np.sqrt(background_sumw2),0),
                                np.append(background_sum+np.sqrt(background_sumw2),0),
                                step='post', color='lightgray', alpha=0.5, hatch='////')            
        if len(data_hists)>0:           
            mh.histplot(data_hists,data_edges[0],
                        histtype='errorbar',
                        stack=False,
                        label=data_labels,
                        color=data_colors,
                        sort='label',
                        w2method='poisson',
                        ax=ax
                        )

        #then stack backgrounds and then draw band
        ax.legend(loc=legend_loc)
        mh.cms.label(self.label, data=True, lumi=self.lumi, ax=ax, loc=0)
        #fix the lower limit
        lims=plt.ylim()
        plt.ylim(0.0, lims[1])
        if xlabel!="":
            ax.set_xlabel(f"{xlabel} ({xunits})")
        if self.stack:
            ax.set_ylabel("Events")
        else:
           ax.set_ylabel("a.u")
             
        if show:
            plt.show()
            

        
        
                    
            
        



        
