clear
clc
load('TransformationMatsForAuCd.mat');
[mHats, bs] = SolveTwinningAndInterfaceEqs(Us);
mHatsUnique = FindUniqueFamilies(mHats);
save('HabitPlanesForAuCd', 'mHatsUnique')
