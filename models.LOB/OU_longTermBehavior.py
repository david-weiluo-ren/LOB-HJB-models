

from abstractLOB.poisson_OU_implicit import Poisson_OU_implicit
from pylab import *
import numpy as np
rcParams['figure.figsize'] = 15, 8

myObj_largeAlpha_largeA_noBeta_longTermMean3_halfLongTime = Poisson_OU_implicit(A=100, s_0=3, kappa=1.5, sigma_s=6.0, num_time_step=5000,\
                                  half_S=5, delta_t=0.001, N=20.0, half_I_S=200,\
                                  beta=-0, q_0=0.0, alpha=40,gamma=1.0, s_long_term_mean=3.0,\
                                   new_weight=0.02, lambda_tilde=0.2,iter_max=2000,\
                                  abs_threshold_power = -6, rlt_threshold_power = -5)
myObj_largeAlpha_largeA_noBeta_longTermMean3_halfLongTime.run()
for _ in xrange(10000):
    tmp = myObj_largeAlpha_largeA_noBeta_longTermMean3_halfLongTime.simulate_forward()
    if tmp[0]:
        break
plot(np.linspace(0, 1, len(myObj_largeAlpha_largeA_noBeta_longTermMean3_halfLongTime.simulate_price_a)), myObj_largeAlpha_largeA_noBeta_longTermMean3_halfLongTime.s[:-1])
plot(np.linspace(0, 1, len(myObj_largeAlpha_largeA_noBeta_longTermMean3_halfLongTime.simulate_price_a)), myObj_largeAlpha_largeA_noBeta_longTermMean3_halfLongTime.simulate_price_a)
plot(np.linspace(0, 1, len(myObj_largeAlpha_largeA_noBeta_longTermMean3_halfLongTime.simulate_price_a)), myObj_largeAlpha_largeA_noBeta_longTermMean3_halfLongTime.simulate_price_b)
xlabel("time")
ylabel("prices")
savefig('longTime_prices.pdf')
plot(np.linspace(0, 1, len(myObj_largeAlpha_largeA_noBeta_longTermMean3_halfLongTime.simulate_price_a)), myObj_largeAlpha_largeA_noBeta_longTermMean3_halfLongTime.q[:-1])
xlabel("time")
ylabel("inventory")
savefig('longTime_inventory.pdf')


