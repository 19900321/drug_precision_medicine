import dataset

# calculate the rate of survive patients at specific time point

def cal_survive_rate(data, ttcos):
    sur_pd = data[data.ttcos <= ttcos]
    sur_rate = sur_pd[sur_pd.censos==1].shape[0]/sur_pd.shape[0]
    return sur_rate


# consider both censors and time to descript the score of each patients,
# the longer time and survive, the score higher
# the shorter time and die is lower
# those short time and survive is expected, score 0, the long time and die is also expected, scroe 0
def cal_survive_score(data, censos, ttcos):
    survive_rate = cal_survive_rate(data, ttcos)
    if censos == 0:
        survive_score = -survive_rate
    elif censos == 1:
        survive_score = 1/survive_rate-1
    return survive_score

def main():
    data = dataset.load_survive_patients()
    data['survive_score'] = data.apply(lambda x: cal_survive_score(data, x.censos, x.ttcos), axis=1)

