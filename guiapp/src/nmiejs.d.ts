export default function nmiejs(): Promise<nmieModule>;


declare interface nmieModule {
    nmie: new () => nmie_class;
}


export class nmie_class {
    constructor(path?: string);
    SetWavelength(wavelength: number): void;
    AddTargetLayerReIm(layer_width: number,
                       re_layer_index: number,
                       im_layer_index: number): void;
    SetModeNmaxAndType(mode_n: number, mode_type: number): void;
    ClearTarget(): void;
    RunMieCalculation(): void;
    RunFieldCalculationPolar(outer_arc_pos: number,
                             radius_pos: number,
                             from_Rho: number,   to_Rho: number,
                             from_Theta: number,   to_Theta: number,
                             from_Phi: number,   to_Phi: number,
                             isIgnoreAvailableNmax: number): void;
    GetFieldEabs(): number[];
    GetQsca(): number;
    GetQabs(): number;
    GetQext(): number;
}

