# SHG-Light-Amplitude-Versus-Phase-Difference-Simulator
A little program aiming at simulating the behavior of the SHG light stimulated from crystals.

The program I wrote outputs the schemes of SHG light amplitude versus phase difference.
Below is the first scheme under different method of phase matching. 
The black dashed line represents the theoretically perfect phase matching condition.
The orange line represents the theoretically condition while using QPM method.
The yellow line represents the theoretically condition while using APP (our intention) method.
The purple line represents the theoretically completely mismatched condition.
<div style="text-align: center;">
  <img src="https://github.com/user-attachments/assets/d84905cc-3a74-4a19-a459-eb315b245ca5" width="381" height="192" alt="image">
</div>

Herein I simulated the curve trend of APP method under different error type and different error rate.
The error type includes independent, that the errors don’t accumulate and dependent, that the errors accumulate. More figurative definition graphs are showed on the right side. 

For the independent mode the center of every domain is fixed and the errors only happen to every border.

For the dependent mode the center of every domain changes due to the error caused by the previous domains.

It is easy to tell the dependent mode is more unstable than the independent one.

<img src="https://github.com/user-attachments/assets/6cd79bf5-2f7a-4e37-bcd7-e78a730168c2" width="298" height="495" alt="Picture1">


The comparation chart is showed below.
The blue line represents the theoretical condition while using APP method.
The orange line represents the condition with a 20% independent fab error while using APP method.
The yellow line represents the condition with a 20% dependent fab error while using APP method.
The purple line represents the condition with a 40% independent fab error while using APP method.
The green line represents the condition with a 40% dependent fab error while using APP method.
We can see that the dependent mode is way more unstable than the independent one. But unfortunately, the fabrication process is always dependent with the former steps in practical. We should try our best to decrease the fabrication error.

<div style="text-align: center;">
  <img width="499" height="300" alt="image" src="https://github.com/user-attachments/assets/653f9ca4-b75a-4ad0-8419-28e33e496207" />
</div>

I also calculated the average final amplitude of different APP methods after the phase difference of 20𝜋 with 100+ runs. The result is showed below. You can find this simulation code in my repo named as SHG_mutiple_average_simulator

<div style="text-align: center;">
  <img width="502" height="292" alt="image" src="https://github.com/user-attachments/assets/0e19dd07-d729-4a08-a8bd-c687bcdae5d1" />
</div>
