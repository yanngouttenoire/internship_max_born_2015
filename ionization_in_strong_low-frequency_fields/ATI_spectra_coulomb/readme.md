# Friday 8/05: I didn't manage to reproduce result in the article

 I can't get the three typical regions:

  - We can briefly see the region I and the expected end at $10 \times U_{p} = 5$ is actually at $2.5 = 5 \times u_{p}$ 

 I tried to :

  - decrease the step dt from 1 to 0.01 -> nothing
  - increase nField from 1000 to 100 0000 and nVPerp from 100 to 2000 -> the curve is smoother and the end at 2.5 is best defined
  - increase the ionization threashold from 100 ua to 600 ua
  - decrease the coulomb threashold below which the coulomb potential is switched off -> the spectra is fully modified
  
Ideas:

 A. Replace the Coulomb Potential $-\frac{1}{r}$ by the Soft Coulomb Potential $-\frac{1}{\sqrt{r^{2}+a}}$
 
  - At nFied, nVPerp and dt fixed, the region of high velocities is smoother

 B. Let's try to plot the trajectories and compare them to the ones in the article

  - splot, 1:2:3, cout<<Vperp<<endl, Vperp=0.123 au (setInitialVPerp), t0=8.050 au (seTBirth)
  - we do not get the same trajectories than in the article
  - for region I, ESS (newton-raphson verified true) gives rhoBirth=15.797 but manually, we get the same behavior than article region I for rhoBirth<12 approximately.
  - Influence of nField: when nField is raised up, the curve is smoother
  - Influence of nVPerp: when nVPerp is increased, the curve seems to morph the one of the article a little better
   - Influence of dt: when we decrease dt, we also have to decrease the soft parameter or else the simulation stops, probably caused to divergence. Elsewhere, as regards the curve transformation, there is not big changes when we decrease dt 

 
  
