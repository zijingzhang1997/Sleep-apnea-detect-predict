---
title: "Sleep apnea detection and prediction based on bed-integrated RF sensor "

---
[Slides](http://zijingzhang1997.github.io/files/sleep/sleep_study_intro.pdf)

<img src='/images/sleep/pic2.png'>  <br/>

**Objective**: Respiratory disturbances during sleep are a prevalent health condition that affects a large adult population. The gold standard to evaluate sleep disorders including apnea is overnight polysomnography, which requires a trained technician for live monitoring and post-processing scoring. Currently, the disorder events can hardly be predicted using the respiratory waveforms preceding the events.  The objective of this work is to develop an autonomous system to detect and predict respiratory events reliably based on real-time covert sensing. 


**Methods**: A bed-integrated radio-frequency (RF) sensor by near-field coherent sensing (NCS) was employed to retrieve continuous respiratory waveforms without user’s awareness. Overnight recordings were collected from 27 patients in the Weill Cornell Center for Sleep Medicine. We extracted respiratory features to feed into the random-forest machine learning model for disorder detection and prediction. The technician annotation, derived from observation by polysomnography, was used as the ground truth during the supervised learning. 

<img src='/images/sleep/pic3.png'>  <br/>


<img src='/images/sleep/pic4.png'>  <br/>
Waveform examples from NCS and PSG in the epochs labelled as (a) normal; (b) OSA; (c) hypopnea. 


**Results**: Apneic event detection achieved a sensitivity and specificity up to 88.6% and 89.0% for k-fold validation, and 83.1% and 91.6% for subject-independent validation.  Prediction of forthcoming apneic events could be made up to 90 s in advance. Apneic event prediction achieved a sensitivity and specificity up to 81.3% and 82.1% for k-fold validation, and 80.5% and 82.4% for subject-independent validation. The most important features for event detection and prediction can be assessed in the learning model.  A bed-integrated RF sensor can covertly and reliably detect and predict apneic events. 


<img src='/images/sleep/pic5.png'> 
<img src='/images/sleep/pic6.png'>  <br/>
Detect and predict sleep disorders using ML model.

**Significance**: Predictive warning of the sleep disorders in advance can intervene serious apnea, especially for infants, servicemen, and patients with chronic conditions. 
  
