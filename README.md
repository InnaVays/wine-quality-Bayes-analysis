# wine-quality-Bayes-analysis

## Project Overview
The goal of this project is to perform a Bayesian analysis on the chemical characteristics of wine and their effects on the wine quality score. Two datasets were examined, containing information about the same chemical parameters for the red and white variants of the Portuguese “Vinho Verde” wine, as well as the wine quality score. The score is on a scale from 0 to 10, with 10 being considered "outstanding". We implemented inference for 3 model. Explored nonlinear transformations and interactions between original explanatory variables, and sampled from posterior predictive distributions. 

Three models were developed and compared in this project:

- Logistic Regression
- Ridge Logistic Regression with variable transformations
- Multinomial regression

## Summary
In this project, we performed a Bayesian analysis of the physicochemical properties of red and white wines using two datasets of the Portuguese “Vinho Verde” wine. 

For each model, we obtained a posterior predictive distribution for a test dataset and estimated its predictability by comparing the results to the observed quality score. However, despite the efforts to tune parameters, the best predictive power was shown by the simplest model with an intercept.

Both cases demonstrated that the intercept takes a value of a dominant category, suppressing all other parameters. This approach proved quite effective for our datasets, where most of the observed values lay in central categories. The poorest result was obtained for multinomial models, as all models had a tendency to choose just one dominant category and misrepresent the others.
