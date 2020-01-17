#pragma once
#include "transition_data.h"

// the functions add the specified transitions to the transition data list:
void ch3oh_classI_trans_list(const energy_diagram*, std::list<transition>& trans_list);
void ch3oh_classII_trans_list(const energy_diagram*, std::list<transition>& trans_list);

void nh3_trans_list(const energy_diagram*, std::list<transition>& trans_list);
void oh_trans_list(const energy_diagram*, std::list<transition>& trans_list);