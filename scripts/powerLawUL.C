// computing the upper limit on the normalisation of a PowerLaw
Double_t powerLaw(Double_t *x, Double_t *par)
{
	// x = log10(E_true / GeV) [gLike histos are in log10(E / GeV)]
	// par[0] = normalisation
	// par[1] = reference
	// par[2] = index (positive)
	Double_t energy = TMath::Power(10, x[0]);
	Double_t e_ratio = energy / par[1];
	Double_t f = par[0] * TMath::Power(e_ratio, - par[2]);
	return f;
}

void powerLawUL()
{
	// initialise the unbinned likelihood only with the input file
	Iact1dUnbinnedLkl* unbn = new Iact1dUnbinnedLkl("inputfile=../data/genericIact_dataIRF_01.root");
	// get the effective time and the effective area's number of bins and energy limits
	TH1F* aeff = (TH1F*) unbn->GetHAeff();
	Double_t obsTime = unbn->GetObsTime();
	Int_t nBinsEtrue = aeff->GetNbinsX();
	Double_t minLog10Etrue = aeff->GetXaxis()->GetXmin();
	Double_t maxLog10Etrue = aeff->GetXaxis()->GetXmax();
	// define the PowerLaw function with the same energy range of the effective area
	TF1 *powerLawTF1 = new TF1("power law function", powerLaw, minLog10Etrue, maxLog10Etrue, 3);
	// set the number of points used to Draw the histogram 
	// using the same number of bins in aeff (to use TH1->Multiply with effective area)
	powerLawTF1->SetNpx(nBinsEtrue); 
	// set the parameters of the PowerLaw
	Double_t phi_0 = 1e-15; // GeV-1 cm-2 s-1
	Double_t E_0 = 300; // GeV
	Double_t index = 3.;
	powerLawTF1->SetParameter(0, phi_0);
	powerLawTF1->SetParameter(1, E_0);
	powerLawTF1->SetParameter(2, index);
	// draw the Power Law
	TCanvas* c1 = new TCanvas("c1", "", 800, 600);
	c1->SetLogy();
	powerLawTF1->Draw();
	powerLawTF1->GetHistogram()->GetYaxis()->SetTitle("d#phi / dE / (GeV^{-1} cm^{-2} s^{-1})");
	powerLawTF1->GetHistogram()->GetXaxis()->SetTitle("log10(E / GeV)");
	powerLawTF1->SetLineWidth(4);
	gPad->Update();
	
	// fetch the histogram from the PowerLaw's TF1
	TH1F* powerLawTH1F = (TH1F*) powerLawTF1->GetHistogram();
	// create a dNdE histogram to feed to gLike
	TH1F* hdNdEinput = new TH1F("hdNdEinput", "dphi / dE multiplied by effective area and time", nBinsEtrue, minLog10Etrue, maxLog10Etrue);
	hdNdEinput->Multiply(aeff, powerLawTH1F, obsTime);
	// draw the dNdE histogram
	TCanvas* c2 = new TCanvas("c2", "", 800, 600);
	c2->SetLogy();
	hdNdEinput->GetYaxis()->SetTitle("dN / dE / (GeV^{-1})");
	hdNdEinput->GetXaxis()->SetTitle("log10(E / GeV)");	
	hdNdEinput->Draw("");

	// set the power law as the dNdE histogram for gLike
	unbn->SetdNdESignal(hdNdEinput);
	unbn->PrintData();

	// take integral of the signal convolved with the IRFs
	Double_t dNdEpSignalIntegral = unbn->GetdNdEpSignalIntegral();
	// set the proper units of g 
	unbn->SetUnitsOfG(1. / (obsTime * dNdEpSignalIntegral));
	// set error and run the likelihood maximisation
	unbn->SetErrorDef(2.71);
	unbn->ComputeLklVsG();
	unbn->PrintOverview();
	//unbn->PlotHistosAndData();
	cout << "-> phi_0 = " << unbn->GetGLklMin() << endl;
	TCanvas* c3 = new TCanvas("c3", "", 800, 600);
	unbn->GetLklVsG()->Draw(); // plot the -2logL vs g curve
}