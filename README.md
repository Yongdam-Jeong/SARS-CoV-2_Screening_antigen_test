# Safe return to schools and workplaces: A modelling study to define guidelines for antigen screening in schools and workplaces to mitigate COVID-19 outbreaks



Yong Dam Jeong, Keisuke Ejima, Kwang Su Kim, Shoya Iwanami, William Hart, Robin N. Thompson, Il Hyo Jung, Shingo Iwami, Marco Ajelli, Kazuyuki Aihara


The attached code can compute simulations for screening using antigen tests in schools and workplaces.
It is flexible to change conditions for the screening simulation in the code.

1) 'populationParameters_sym.txt' and 'populationParameters_asym.txt' are the estimated parameters of viral dynamics model for symptomatic and asymptomatic SARS-CoV-2 patients, respectively.
2) 'Baseline.R' is a code for simulating screening schedules at a facility under a baseline epidemiological setting (i.e., reproduction number = 2 & asymptomatic ratio = 70%).
3) 'Sensitivity_epidemic.R' is a code for simulating screening schedules at a facility under various epidemiological settings (i.e., reproduction number & asymptomatic ratio).
4) 'Sensitivity_resources.R' is a code for simulating Schedule 1 at a facility under various epidemiological settings (i.e., reproduction number & asymptomatic ratio) and resources (i.e., sensitivity of antigen test, screening initiation timing, and available tests per person).

* Text files and R codes above should be in the same location.
