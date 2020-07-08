function [hydroMaxError katieHMaxError fullMaxError] ...
    = getErrorsFromApproxData(expNum,h,etaFilter,etaHydro,etaKatieH,etaFull)

etaFilterMax = max(etaFilter);
etaHydroMax = max(etaHydro); hydroMaxError = abs(etaFilterMax - etaHydroMax)/etaFilterMax/.01;
etaKatieHMax = max(etaKatieH); katieHMaxError = abs(etaFilterMax - etaKatieHMax)/etaFilterMax/.01;
etaFullMax = max(etaFull); fullMaxError = abs(etaFilterMax - etaFullMax)/etaFilterMax/.01;
