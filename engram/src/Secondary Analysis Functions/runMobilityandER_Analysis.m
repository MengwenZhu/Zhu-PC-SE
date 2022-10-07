function Data_Summary = runMobilityandER_Analysis(Inscopix1, Inscopix2, Noldus1, Noldus2, Start_T1, Start_T2, PlaceCell_Summary1, PlaceCell_Summary2, dateToFind, animalName, drug, session1_dose, session2_dose, session1_context, session2_context, genotype, exptParadigm, ds, ap)
% runEngram_Cell_Analysis - calls Engram_Cell_Analysis_V1 and appends fields

Data_Summary = Mobility_and_ER_Analysis(Inscopix1, Inscopix2, Noldus1, Noldus2, Start_T1, Start_T2, PlaceCell_Summary1, PlaceCell_Summary2, ds, ap);
% append fields to structure Data_Summary created above
Data_Summary.dateToFind = dateToFind;
Data_Summary.animalName = animalName;
Data_Summary.drug = drug;
Data_Summary.session1_dose = session1_dose;
Data_Summary.session2_dose = session2_dose;
Data_Summary.session1_context = session1_context;
Data_Summary.session2_context = session2_context;
Data_Summary.genotype = genotype;
Data_Summary.exptParadigm = exptParadigm;

end