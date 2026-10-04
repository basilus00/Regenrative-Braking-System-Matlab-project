function [Iq_ref,Id_ref] = get_mtpa_vals(Id_ref_mtpa,Iq_ref_mtpa)
    Is_ref = sqrt(Id_ref_mtpa^2 + Iq_ref_mtpa^2);
    if Is_ref > I_max
        Id_ref_mtpa = Id_ref_mtpa * I_max / Is_ref; %assure the references are always within the limits (limitation)
        Iq_ref_mtpa = Iq_ref_mtpa * I_max / Is_ref;
    end
    Id_ref = Id_ref_mtpa;
    Iq_ref = Iq_ref_mtpa;
end
