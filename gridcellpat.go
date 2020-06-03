package main

import (
	"math"
	"math/rand"
	"os"
	"strconv"

	gridcell "github.com/juishiang/GridCell"

	"github.com/cpmech/gosl/plt"
	"github.com/emer/etable/etable"
	"github.com/emer/etable/etensor"
	"github.com/goki/gi/gi"
)

var test gridcell.Grid_layer
var data [][][]float64

func main() {
	placedevsizeS := []string{"032", "030", "028", "026"} //{"030", "025", "020", "015", "010"}
	placedevsizeF := []float64{0.32, 0.3, 0.28, 0.26}     //{0.3, 0.25, 0.2, 0.15, 0.1}
	norm := []bool{false, true}
	normS := []string{"Nnorm", "Norm"}
	SpRan := []bool{false, true}      //
	SpRanS := []string{"std", "Rand"} //
	var spacesize float64 = 10
	var neunum int64 = 125
	stepsize := 0.1
	totalstep := int((spacesize / stepsize) + 1)
	data = make([][][]float64, totalstep)
	for i := range data {
		data[i] = make([][]float64, totalstep)
		for j := range data[i] {
			data[i][j] = make([]float64, neunum)
		}
	}
	//////initialize
	for idx, placedevsizeFvalue := range placedevsizeF {
		for Ridx, tr := range SpRan {
			for idxN, Nor := range norm {
				test.Init(neunum)
				str, _ := os.Getwd()
				str = str + "/fi_0601/" + placedevsizeS[idx] + SpRanS[Ridx] + normS[idxN]
				csvname := gi.FileName(placedevsizeS[idx] + SpRanS[Ridx] + normS[idxN] + ".csv")
				for i := range test.Grlay {
					if tr {
						test.Grlay[i].Init(i%5+3.0, rand.Float64(), rand.Float64(), spacesize, 0.8, 0.6, placedevsizeFvalue)
					} else {
						test.Grlay[i].Init(i%5+3.0, (float64((i/5)%5) / 5.0), (float64(i/25) / 5.0), spacesize, 0.8, 0.6, placedevsizeFvalue)
					}
					//fmt.Println(rand.Float64(), rand.Float64())
					//fmt.Println(float64(i%5)/5.0, float64(i/25)/5.0)
				}
				for k := range test.Grlay {
					for i := 0; i < totalstep; i++ {
						for j := 0; j < totalstep; j++ {
							data[i][j][k] = test.Grlay[k].Fireact(float64(i)*stepsize, float64(j)*stepsize)
						}
					}
				}
				///////set position array using in the figure//////
				var X [][]float64
				var Y [][]float64
				var Neurondata [][]float64
				X = make([][]float64, len(data))
				Y = make([][]float64, len(data))
				Neurondata = make([][]float64, len(data))
				for i := 0; i < len(data); i++ {
					X[i] = make([]float64, len(data[0]))
					Y[i] = make([]float64, len(data[0]))
					Neurondata[i] = make([]float64, len(data[0]))
					for j := 0; j < len(data[0]); j++ {
						X[i][j] = stepsize * float64(i)
						Y[i][j] = stepsize * float64(j)
					}
				}
				////////end position
				//Normalize the input pattern for each site
				if Nor {
					for i := 0; i < len(data); i++ {
						for j := 0; j < len(data[0]); j++ {
							totalsum := float64(0)
							for k := 0; k < len(data[0][0]); k++ {
								totalsum += data[i][j][k] * data[i][j][k]
							}
							for k := 0; k < len(data[0][0]); k++ {
								data[i][j][k] = math.Sqrt((data[i][j][k] * data[i][j][k]) / totalsum)
								data[i][j][k] *= float64(6)
							}
						}
					}
				}

				////////put the data to etable and tensor
				dt := etable.NewTable("data_test")
				dt.AddRows(totalstep * totalstep)
				//Tdata := make([][]float64, neunum)
				Tetsr1 := etensor.NewStringShape(etensor.NewShape([]int{int(totalstep * totalstep)}, nil, nil))
				Tetsr2 := etensor.NewStringShape(etensor.NewShape([]int{int(totalstep * totalstep)}, nil, nil))
				for i := 0; i < totalstep; i++ {
					for j := 0; j < totalstep; j++ {
						Tetsr1.Set([]int{i*totalstep + j}, "_D:")
						Tetsr2.Set([]int{i*totalstep + j}, "place: ("+strconv.FormatFloat(stepsize*float64(i), 'E', -1, 64)+","+strconv.FormatFloat(stepsize*float64(j), 'E', -1, 64)+")")
					}
				}
				dt.AddCol(Tetsr1, "_H:")
				dt.AddCol(Tetsr2, "Name")
				for k := range test.Grlay {
					etsr := etensor.NewFloat32([]int{int(totalstep * totalstep)}, nil, nil)
					for i := 0; i < totalstep; i++ {
						for j := 0; j < totalstep; j++ {
							//Tdata[k][i*totalstep+j] = data[i][j][k]
							etsr.Set([]int{i*totalstep + j}, float32(data[i][j][k]))
						}
					}
					//etsr := etensor.NewFloat64Shape(etensor.NewShape([]int{int(totalstep * totalstep)}, nil, nil), Tdata[k])
					dt.AddCol(etsr, "Input[2:"+strconv.Itoa(k/25)+","+strconv.Itoa(k%25)+"]")
				}
				dt.ColNames[2] = dt.ColNames[2] + "<2:25,5>"
				dt.SaveCSV(csvname, 0, true)

				//////////
				for k := 0; k < len(data[0][0]); k++ {
					for i := 0; i < len(data); i++ {
						for j := 0; j < len(data[0]); j++ {
							Neurondata[i][j] += data[i][j][k]
						}
					}
				}
				plt.Reset(false, nil)
				levels := []float64{0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1}
				plt.ContourF(X, Y, Neurondata, &plt.A{CmapIdx: 4, Nlevels: 15 /*, Levels: levels*/})
				plt.Gll("x", "y", nil)
				plt.Equal()
				filename := "plt_contour_all" // + strconv.Itoa(k)
				plt.Save(str, filename)

				for k := 0; k < len(data[0][0]); k++ {
					for i := 0; i < len(data); i++ {
						for j := 0; j < len(data[0]); j++ {
							Neurondata[i][j] = data[i][j][k]
						}
					}
					plt.Reset(false, nil)
					plt.ContourF(X, Y, Neurondata, &plt.A{CmapIdx: 4, Levels: levels})
					plt.Gll("x", "y", nil)
					plt.Equal()
					filename := "plt_contour" + strconv.Itoa(k)
					plt.Save(str, filename)
				}
				levels = []float64{0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4}
				for i := 0; i < len(data); {
					for j := 0; j < len(data[0]); {
						for a := 0; a < len(data); a++ {
							for b := 0; b < len(data[0]); b++ {
								ip := float64(0)
								for k := 0; k < len(data[0][0]); k++ {
									ip += data[a][b][k] * data[i][j][k]
								}
								Neurondata[a][b] = ip
							}
						}
						plt.Reset(false, nil)
						plt.ContourF(X, Y, Neurondata, &plt.A{CmapIdx: 4, Nlevels: 20 /*Levels: levels*/})
						plt.Gll("x", "y", nil)
						plt.Equal()
						filename := "ip_" + strconv.Itoa(i) + "_" + strconv.Itoa(j)
						plt.Save(str+"/ip", filename)
						j += 5
					}
					i += 5
				}
			}
		}
	}
}
