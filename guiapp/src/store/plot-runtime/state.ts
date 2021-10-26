export interface spectraData {
  WLs: number[]
  Qsca:number[]
  Qabs:number[]
  Qext:number[]
  Qsca_n:number[][][]
  Qabs_n:number[][][]
  Qext_n:number[][][]
}

export interface plotRuntimeStateInterface {
  WLs: number[]
  Qsca:number[]
  Qabs:number[]
  Qext:number[]
  Qsca_n:number[][][]
  Qabs_n:number[][][]
  Qext_n:number[][][]
}

function state(): plotRuntimeStateInterface {
  const WLs:number[] = []
  const Qsca:number[] = [], Qabs:number[] = [], Qext:number[] = []
  const Qsca_n:number[][][] = [[], []]
  const Qabs_n:number[][][] = [[], []]
  const Qext_n:number[][][] = [[], []]
  return { WLs,
    Qsca, Qabs, Qext, Qsca_n, Qabs_n, Qext_n
  }
}

export default state
