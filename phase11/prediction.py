from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.metrics import r2_score,mean_absolute_error,mean_squared_error
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results'
def main():
    policy=json.loads((P/'prediction_policy.json').read_text());features=policy['predictors'];x=pd.read_csv(R/'new_descriptor_contrasts.csv');x=x[x.metric.isin(features)].pivot(index=['gene','batch'],columns='metric',values='effect').reset_index();y=pd.read_csv(ROOT/'phase6/results/adaptive_conditional_contrasts.csv')[['gene','batch','effect']].rename(columns={'effect':'alignment'});data=x.merge(y,on=['gene','batch'],validate='one_to_one');assert len(data)==54 and data[features+['alignment']].notna().all().all();predictions=[];foldrows=[];summaries=[]
    for grouping in ['gene','batch']:
        groups=data[grouping].to_numpy();X=data[features].to_numpy();Y=data.alignment.to_numpy();pred=np.empty(len(data))
        for train,test in LeaveOneGroupOut().split(X,Y,groups):
            losses=[]
            for alpha in [.1,1.,10.,100.]:
                err=[]
                for ti,vi in LeaveOneGroupOut().split(X[train],Y[train],groups[train]):
                    model=make_pipeline(StandardScaler(),Ridge(alpha=alpha));model.fit(X[train][ti],Y[train][ti]);err.append(mean_squared_error(Y[train][vi],model.predict(X[train][vi])))
                losses.append(float(np.mean(err)))
            alpha=[.1,1.,10.,100.][int(np.argmin(losses))];model=make_pipeline(StandardScaler(),Ridge(alpha=alpha));model.fit(X[train],Y[train]);pred[test]=model.predict(X[test]);foldrows.append({'outer_grouping':grouping,'held_out':str(groups[test][0]),'train_observations':len(train),'test_observations':len(test),'selected_alpha':alpha,'inner_mean_squared_error':min(losses)})
        q=data.copy();q['outer_grouping']=grouping;q['predicted_alignment']=pred;predictions.append(q);summaries.append({'outer_grouping':grouping,'observations':len(Y),'groups':len(set(groups)),'out_of_fold_R2':float(r2_score(Y,pred)),'out_of_fold_MAE':float(mean_absolute_error(Y,pred))})
    pd.concat(predictions,ignore_index=True).to_csv(R/'shape_to_alignment_OOF_predictions.csv',index=False);pd.DataFrame(foldrows).to_csv(R/'shape_to_alignment_folds.csv',index=False);pd.DataFrame(summaries).to_csv(R/'shape_to_alignment_prediction_summary.csv',index=False);print(json.dumps(summaries,indent=2))
if __name__=='__main__':main()
