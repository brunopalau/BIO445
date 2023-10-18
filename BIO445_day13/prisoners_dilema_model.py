#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 11:25:20 2022

@author: bp
"""

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import SGDClassifier
from sklearn import svm
import random



    


def payout(action_p1, action_p2, pay = np.array([5, 3, 1, 0])):

    if action_p1 and action_p2:
        return np.array([pay[1], pay[1]])

    if action_p1 and not action_p2:
        return np.array([pay[3], pay[0]])

    if not action_p1 and action_p2:
        return np.array([pay[0], pay[3]])

    if not action_p1 and not action_p2:
        return np.array([pay[2], pay[2]])


def add_feature(X, partnerpastactions, ownpastactions):
    if len(ownpastactions) < 10:
        n = len(ownpastactions)
        start = 1
        
        # create feature for this round
        # correlation = np.corrcoef(np.vstack((ownpastactions,partnerpastactions)))[0,1]
        total_true = np.sum(partnerpastactions[start:])
                
    else:
        start = len(partnerpastactions)-10 
        # create feature for this round
        # correlation = np.corrcoef(np.vstack((ownpastactions,partnerpastactions)))[0,1]
        total_true = np.sum(partnerpastactions[start:])
        
    after_true = np.sum(partnerpastactions[list(np.argwhere(ownpastactions[0:-2] == True).ravel()+1)])
    after_false = np.sum(partnerpastactions[list(np.argwhere(ownpastactions[0:-2] == False).ravel()+1)])
    last_1 = partnerpastactions[-1]
    new_entry = {"Total_true": total_true, "After_false":after_false, "After_true":after_true, "Last_1": last_1}#, "Last_2": last_2, "Last_3": last_3}
    new_entry = pd.DataFrame(data=new_entry,index=[0])
    
    # append new entry
    X = pd.concat([X, new_entry])

    return X
    


def simulate_2P_PD(strategy_p1, strategy_p2, N, model1, model2, pay = np.array([5, 3, 1, 0])):
    # create dataframe to store values
    d = { "Total_true": [], "After_false": [], "After_true": [], "Last_1": []}#, "Last_2": [], "Last_3": []}
    X_for_1 = pd.DataFrame(data = d)
    X_for_2 = pd.DataFrame(data = d)
    
    predictions_for_2 = np.array([])
    predictions_for_1 = np.array([])

    total_pay = np.array([0, 0])

    past_actions_p1 = np.array([])
    past_actions_p2 = np.array([])
    
    

    for i in range(N):
        prediction_for_1 = np.NaN
        prediction_for_2 = np.NaN
        
        
        if i > 50:
            if not all(past_actions_p1 == True) and not all(past_actions_p2 == False):
                model1.fit(X_for_1[:-1],past_actions_p1[1:])
                prediction_for_1 = model1.predict(X_for_1.tail(1))
            
            if not all(past_actions_p2 == True) and not all(past_actions_p2 == False) :
                model2.fit(X_for_2[:-1],past_actions_p2[1:])
                prediction_for_2 = model2.predict(X_for_2.tail(1))
            

        out1 = strategy_p1(past_actions_p1, past_actions_p2)
        out2 = strategy_p2(past_actions_p2, past_actions_p1)

        total_pay += payout(out1, out2)

        past_actions_p1 = np.append(past_actions_p1, out1)
        past_actions_p2 = np.append(past_actions_p2, out2)
        
        # create features and save
        X_for_1 = add_feature(X_for_1, past_actions_p1, past_actions_p2)
        X_for_2 = add_feature(X_for_2, past_actions_p2, past_actions_p1)
        
        # store predictions
        predictions_for_2 = np.append(predictions_for_2,prediction_for_2)
        predictions_for_1 = np.append(predictions_for_1,prediction_for_1)
        
        
    # compute accuracy
    accuracy_for_1 = np.sum(predictions_for_1 == past_actions_p1)/(len(predictions_for_1) - sum(np.isnan(predictions_for_1)))
    accuracy_for_2 = np.sum(predictions_for_2 == past_actions_p2)/(len(predictions_for_2) - sum(np.isnan(predictions_for_2)))
    

    return total_pay, accuracy_for_1, accuracy_for_2

def tournament(list_of_strategies, N, Nsim, model1, model2, pay = np.array([5, 3, 1, 0])):
    n = len(list_of_strategies)
    total_payouts = np.zeros((n ,n))
    av_payout_mat = np.zeros((n ,n))
    av_accuracy_mat = np.zeros((n ,n))

    for i in range(n):
        for j in range(i+1, n):
            sum_payout = np.array([0, 0])
            sum_accuracy_1 = 0
            sum_accuracy_2 = 0
            for sim in range(Nsim):
                total_pay, accuracy_for_1, accuracy_for_2 = simulate_2P_PD(list_of_strategies[i], list_of_strategies[j], N, model1, model2)
                sum_accuracy_1 += accuracy_for_1
                sum_accuracy_2 += accuracy_for_2
                sum_payout += total_pay
                
            # avg performance per battle
            average_payout = sum_payout / Nsim
            total_payouts[i] += average_payout[0]
            total_payouts[j] += average_payout[1]
            
            # avg accuracy
            average_accuracy_1 = sum_accuracy_1/Nsim
            average_accuracy_2 = sum_accuracy_2/Nsim
            av_accuracy_mat[i, j] = average_accuracy_1
            av_accuracy_mat[j, i] = average_accuracy_2

            # append to matrix
            av_payout_mat[i, j] = average_payout[0]
            av_payout_mat[j, i] = average_payout[1]


    return total_payouts, av_payout_mat, av_accuracy_mat



def tit_for_tat(ownpast, partnerpast):

    n = len(ownpast)

    if n == 0:
        return True

    else:
        return partnerpast[-1]


def tit_for_two_tats(ownpast, partnerpast):

    n = len(ownpast)

    if n == 0:
        return True

    if not partnerpast[-2:].any():
        return False

    else:
        return True

def two_tits_for_tat(ownpast, partnerpast):

    n = len(ownpast)

    if n == 0:
        return True

    if partnerpast[-2:].all():
        return False

    else:
        return True
    
    
    

d = { "Total_true": [], "After_false": [], "After_true": [], "Last_1": []}#, "Last_2": [], "Last_3": []}
features = pd.DataFrame(data = d)
    
total_pay, accuracy_for_1, accuracy_for_2 = simulate_2P_PD(two_tits_for_tat, tit_for_tat, 1000, LogisticRegression(),LogisticRegression())
    
total_payouts, av_payout_mat, av_accuracy_mat = tournament([two_tits_for_tat, tit_for_tat], 500, 10, LogisticRegression(),LogisticRegression())
    

