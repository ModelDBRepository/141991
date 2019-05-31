This is the model for the paper:

Jercog PE, Svirskis G, Kotak VC, Sanes DH, Rinzel J (2010) Asymmetric
excitatory synaptic dynamics underlie interaural time difference
processing in the auditory system. PLoS Biol 8:e1000406

The model code was contributed by Pablo Jercog.

Two Models for synaptic inputs statistics for the MSO neuron.

These two models generate the conductance, membrane potential and low
threshold potassium current time courses for bilateral inputs at
different ITDs. These codes also produce ITD-tuning curves for the
cell model given these synaptic input models.

This Medial Superior Olivary (MSO) neuron model generates responses to
two synaptic input trains that mimic interaural time delay (ITD)
inputs from the two ears for pure tone stimuli. An MSO neuron
generates spikes if the 'bilateral' inputs are nearly coincident in
time, performing a crucial neuronal computation in the mammalian
brainstem for sound localization. The model is a point neuron model
with ionic currents from Rothman & Mannis (2003) and with an update of
the low threshold potassium current (IKLT) measured in-vitro by
Mathews & Jercog et al (2010). This model in conjunction with the
synaptic input models presented here has been used to gain insight
into mechanisms that account for experimentally observed asymmetries
in ITD tuning (Brand et al, 2002). Asymmetry and displacement of the
ITD response function is achieved in the model by the interplay
between asymmetry of the excitatory inputs arriving from the two sides
and the precise voltage dependent activation of IKLT. In Jercog et al
(2010) we propose two different mathematical ways, physiologically
plausible scenarios, of generating the asymmetry in the bilateral
synaptic input events. Here, we present two models for simulating the
stochastic synaptic input trains.

Model 1: The "Gaussian-jittered poisson-modulated synaptic inputs
model" is a simple model that we used in figure 3 (Jercog et al 2010)
for presynaptic event times. It involves a partially
rectified-sinusoidal time-modulated poisson function with additional
gaussian jitter to further randomize each individual event time on a
given cycle. This model best resembles our recorded data (in vitro,
bilaterally stimulated MSO neuron) in terms of variability in
composite EPSP amplitude, half-width and rising slope. A potential
source of gaussian statistical distribution that we include in the
time of presynaptic events is attributed to a possible distribution of
delays of pre-synaptic pathways' propagation-time (say, gaussian
distribution of axonal lengths) and location of pre-synaptic terminals
over the MSO dendrites (also gaussian as the simplest biological
plausible distribution).

Model 2: The "Carney-like" synaptic input model generates presynaptic
input event times using Carney's 1993 auditory nerve model, modifying
the parameter values in order to obtain summated EPSPs with a rising
slope in the range of values as observed in our experiments (figure 2
in Jercog et al 2010). This model was used as a proof of principles to
generate figures 4 and 5 (Jercog et al 2010) based on a model, well
know in the auditory community. See also: Carney 1993, Brughera et al
1996, Rothman et al 1993, first proposed by Johnson et al 1986.

How to run the models:

# gcc [code-name.c] -lm -O3 -o [executable-name]
# ./[executable-name]

####################################################################

Outputs from these codes:
(Plotting the outputs in gnuplot or matlab it's easy)

"Model1_STA_epsg_epsp_iklt.dat" different averaged values triggered by
the spike detection
first column = ITD
second column = time before the spike's peak (ms)
third column = spike triggered voltage average
forth column = spike triggered voltage lowest-range trace
fifth column = spike triggered voltage highest-range trace
sixth column = non-spike triggered voltage average
seventh column = non-spike triggered voltage lowest-range trace
eighth column = non-spike triggered voltage highest-range trace
ninth column = non-spike triggered voltage std.dev. 
tenth column = spike triggered low-treshold potassium current (IKLT)
               average
eleventh column = spike triggered low-treshold potassium conductance
                  (gKLT) average
twelfth column = non-spike triggered low-treshold potassium
                 conductance (gKLT) average
thirteenth column = non-spike triggered ipsilateral excitatory
                    conductance average
fourteenth column = non-spike triggered contralateral excitatory
                    conductance average
fifteenth column = spike triggered ipsilateral excitatory conductance
                   average
sixteenth column = spike triggered contralateral excitatory
                   conductance average

"Model1_train_EPSGs-EPSPs.dat" is the train of EPSPs and EPSGs for
ITD = 0msec.

first column=time
second column=voltage time course
third column = sinusoidal time modulation of the synaptic poisson rate
forth column = total synaptic input
fifth column = ipsilateral excitatory synaptic input
sixth column = contralateral excitatory synaptic input

"Model1_ITD-curve.dat" ITD response function
first column = ITD
second column = spike count/(total number of sinusoidal cycles that
                modulate the poisson rate)

####################################################################

"Model2_syn-events-ipsi.dat", "Model2_syn-events-contra.dat" are the
output files that contain:
first column=current time of the simulation;
second column= index of the synaptic input that is activated
third column= the exact time where the i-th event is located
forth column= value of the poisson event rate function at the time of
              the decision of creating a synaptic event.
(ipsi and contra correspond to each bilateral input.)

"Model2_train_EPSGs-EPSPs.dat" is the train of EPSPs and EPSGs for
ITD = 0msec.
first column=time
second column=voltage time course
third column = sinusoidal time modulation of the synaptic poisson rate
forth column = total synaptic input
fifth column = ipsilateral excitatory synaptic input
sixth column = contralateral excitatory synaptic input

"Model2_ITD-curve.dat" ITD response function
first column = ITD
second column = spike count/(total number of sinusoidal cycles that
                modulate the poisson rate )

"Model2_STA_epsg_epsp_iklt.dat" different averaged values triggered by
the spike detection
first column = time previous to the spike's peak
second column=spike triggered voltage average previous to the spike's
              peak
third column = spike triggered ipsilateral EPSG average previous to
               the spike's peak
forth column = spike triggered contralateral EPSG average previous to
               the spike's peak
fifth column = spike triggered contralateral IPSG average previous to
               the spike's peak
