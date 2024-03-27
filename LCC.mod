TITLE Big conductance K channel in ER

NEURON {
	SUFFIX LCC
	USEION k READ ki, ko WRITE ik
        RANGE gbar, ki, ik				    
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)			
	FARADAY = 96485 (coul)
	R = 8.3134 (joule/degC)
		
}



PARAMETER {          : parameters that can be entered when function is called in cell-setup
        v               (mV)
        celsius = 34	(degC)
	dt              (ms)
        gbar = 0     (mho/cm2) : initialized conductance
	ki = 6 (mM)
	ko = 120 (mM)				
        }


ASSIGNED {
	ik (mA/cm2)
        
}


BREAKPOINT {
	ik = gbar*ghk(v, ki, ko)

}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f
        f = KTF(celsius)
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {
     KTF = (1e+3)*(R*(celsius + 273.15)/FARADAY/2)
     :KTF = RT/(zF)
}


FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}






















