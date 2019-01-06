simulated:
	sbt ";project examples; runMain examples.SimulateGp; runMain examples.FitGp; runMain examples.ParametersSimulatedGp; runMain examples.PosteriorPredictive"

fit_temp_dlm:
	sbt ";project examples; runMain examples.FitTemperatureDlm"

forecast_temperature_dlm:
	sbt ";project examples; runMain examples.ForecastTemperatureDlm; runMain examples.StateTemperatureDlm; runMain examples.ForecastTestDlm"

fit_temp_gp:
	sbt ";project examples; runMain examples.FitTemperatureResiduals; runMain examples.ForecastGp"

