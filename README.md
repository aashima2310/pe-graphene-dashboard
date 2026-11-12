# PE / Graphene MD Explorer

Extract the ZIP, open a terminal inside pe_dashboard, and run:

    python -m pip install -r requirements.txt
    python -m streamlit run app.py

Use Python 3.11 or 3.12. Keep the data and .streamlit folders beside app.py.

Replace data/results_big.csv or data/results_xyz.csv with updated files using the same columns, then click Reload data. Charts and counts recalculate. Only completed configurations should be in the CSV. Do not add duplicate runs or treat directional tests as separate seeds.

This release contains ONLY the individually uploaded November 12 data: 32 small-system rows, 12 large-system rows, 15 adhesion rows, 48 interface-profile rows, two representative large-system y curves, and the provided OVITO image. The older ZIP is not used.

95% confidence intervals are Student-t intervals of composition means. Bands on interface profiles are supplied standard deviations. No fixed p-values, claimed significance, or detection limits are hard-coded. Stress smoothing is display-only. Theory inputs start unset; enter sourced parameters to plot a Halpin–Tsai comparison. The app does not infer a scientifically established optimum from the highest mean.

Theoretical model parameters and raw cooling curves are not supplied. To add a temperature-density curve page later, provide tg.txt for a representative neat and filled system with temperature/density column definitions. Existing Tg estimates are available now.

Team roles and protocol descriptions follow the proposal/project description and should be confirmed before presenting. The dashboard does not independently validate the MD simulation protocol.
