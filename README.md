## Pephemeral
Pephemeral is a project I finished while working as a Health Data Science Fellow at Insight Data Science.  
The purpose of the project is to predict peptide drug stability from user input sequence information. 
The structure of the code is in accordance with the recommended layout for Pyton Flask web app development.  

The web app itself can be found at *pephemeral.com*  

The main code can be found in the following two files:  

#### compute.py:   
trains the machine-learning model from training set, generates feature vector from user sequence input and specified conditions, handles degenerate sequence specification, and finally returns the predicted peptide half-life range in minutes  

#### view.py: 
renders the html web page for user input and model output

