export const flexRowTitleStyle = 'width:10em; margin: auto;'
export const basicSelectorWidthStyle = '6.5em'

export const inputWithUnitsHistoryLength = 5
export const inputWithUnitsTooltipDigits = 2
export const inputWithUnitsInlineDigits = 1

export const inputWithUnitsTitleWidthStyle = '4em'
export const inputWithUnitsBodyWidthStyle = '10em'
export const inputWithUnitsUnitsWidthStyle = '3em'
function getInt(val:string) {return parseInt(val.slice(0,-2))}
export const basicWidthStyle = (
    getInt(inputWithUnitsBodyWidthStyle)
    + getInt(inputWithUnitsTitleWidthStyle)
    + getInt(inputWithUnitsUnitsWidthStyle)
    + 0.71 // To get the same width in GUI as all three parts above joined in inputWithUnits component
).toString()+'em'


export const maxNumberOfModesToPlot = 10
export const maxNumberOfLayers = 10

