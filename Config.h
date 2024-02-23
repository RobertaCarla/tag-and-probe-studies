#ifndef CONFIG_H
#define CONFIG_H

constexpr double kMaxCounts = 1.e5;
constexpr float kMassWindow = 3.f;
constexpr float kSBandsWindow = 5.f;

constexpr int kNbinsx= 24;
constexpr double k_xbins[kNbinsx+1]={0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2., 2.4, 2.8, 3.2, 3.6, 4., 4.4, 5.};  //l'ultimo bin è di 0.6 affinché fosse compatibile con binning precedente
constexpr float epsilon=0.05;


#endif // CONFIG_H