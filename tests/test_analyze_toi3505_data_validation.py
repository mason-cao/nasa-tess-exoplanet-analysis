from __future__ import annotations

import json
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from analyze_toi3505_data_validation import (
    decode_spoc_sector_vector,
    dvt_matches_xml,
    parse_dr122_target_info,
    parse_dv_xml,
    parse_dvt_header,
    parse_exomast_info,
    parse_exomast_table,
)


SYNTHETIC_XML = """<dv:dvTargetResults xmlns:dv="http://www.nasa.gov/2018/TESS/DV"
    ticId="390988385" toiId="3505">
  <dv:planetResults>
    <dv:allTransitsFit fullConvergence="true" modelFitSnr="20.0"
        transitModelName="test_transit_model">
      <dv:modelParameters>
        <dv:modelParameter name="orbitalPeriodDays" value="2.915" uncertainty="0.001"/>
        <dv:modelParameter name="transitEpochBtjd" value="3000.5" uncertainty="0.002"/>
        <dv:modelParameter name="transitDurationHours" value="2.5" uncertainty="0.1"/>
        <dv:modelParameter name="transitDepthPpm" value="3200" uncertainty="200"/>
        <dv:modelParameter name="ratioPlanetRadiusToStarRadius" value="0.06" uncertainty="0.003"/>
        <dv:modelParameter name="minImpactParameter" value="0.8" uncertainty="0.1"/>
      </dv:modelParameters>
    </dv:allTransitsFit>
    <dv:binaryDiscriminationResults>
      <dv:oddEvenTransitDepthComparisonStatistic value="0.25" significance="0.80"/>
    </dv:binaryDiscriminationResults>
    <dv:bootstrapResults significance="1e-20"/>
    <dv:centroidResults>
      <dv:differenceImageMotionResults>
        <dv:msTicCentroidOffsets>
          <dv:meanSkyOffset value="2.0" uncertainty="2.5"/>
        </dv:msTicCentroidOffsets>
        <dv:msControlCentroidOffsets>
          <dv:meanSkyOffset value="3.0" uncertainty="2.5"/>
        </dv:msControlCentroidOffsets>
        <dv:summaryQualityMetric fractionOfGoodMetrics="1.0"/>
      </dv:differenceImageMotionResults>
    </dv:centroidResults>
    <dv:differenceImageResults sector="54" numberOfTransits="8">
      <dv:qualityMetric value="0.95"/>
    </dv:differenceImageResults>
    <dv:ghostDiagnosticResults>
      <dv:coreApertureCorrelationStatistic value="12.0"/>
      <dv:haloApertureCorrelationStatistic value="2.0"/>
    </dv:ghostDiagnosticResults>
    <dv:planetCandidate maxMultipleEventSigma="18.0"
        observedTransitCount="8" expectedTransitCount="9"
        suspectedEclipsingBinary="false">
      <dv:weakSecondary maxMes="2.2" maxMesPhaseInDays="1.1">
        <dv:depthPpm value="300" uncertainty="120"/>
      </dv:weakSecondary>
    </dv:planetCandidate>
  </dv:planetResults>
</dv:dvTargetResults>
"""


class DataValidationTests(unittest.TestCase):
    def test_xml_parser_extracts_selected_report_values(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report.xml"
            path.write_text(textwrap.dedent(SYNTHETIC_XML), encoding="utf-8")
            row = parse_dv_xml(path)

        self.assertEqual(row["sector"], 54)
        self.assertTrue(row["fit_full_convergence"])
        self.assertAlmostEqual(row["fit_depth_ppt"], 3.2)
        self.assertAlmostEqual(row["fit_depth_error_ppt"], 0.2)
        self.assertAlmostEqual(row["odd_even_difference_sigma"], 0.5)
        self.assertAlmostEqual(row["odd_even_planet_favoring_percent"], 80.0)
        self.assertAlmostEqual(row["weak_secondary_depth_ppt"], 0.3)
        self.assertAlmostEqual(row["centroid_to_tic_sigma"], 0.8)
        self.assertAlmostEqual(row["ghost_core_halo_ratio"], 6.0)
        self.assertFalse(row["suspected_eclipsing_binary"])

    def test_dvt_header_parser_reads_summary_keywords(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report_dvt.fits"
            primary = fits.PrimaryHDU()
            primary.header["SECTOR"] = 54
            primary.header["TICID"] = 390988385
            tce = fits.BinTableHDU.from_columns([], name="TCE_1")
            tce.header["TPERIOD"] = 2.915
            tce.header["TEPOCH"] = 3000.5
            tce.header["TDUR"] = 2.5
            tce.header["TDEPTH"] = 3200.0
            tce.header["MAXMES"] = 18.0
            fits.HDUList([primary, tce]).writeto(path)
            row = parse_dvt_header(path)

        self.assertEqual(row["sector"], 54)
        self.assertEqual(row["tic_id"], 390988385)
        self.assertAlmostEqual(row["dvt_depth_ppt"], 3.2)
        self.assertAlmostEqual(row["dvt_maximum_multiple_event_statistic"], 18.0)

    def test_xml_and_dvt_consistency_check_has_clear_tolerances(self) -> None:
        xml_row = {
            "sector": 54,
            "tic_id": 390988385,
            "fit_period_days": 2.915,
            "fit_epoch_btjd": 3000.5,
            "fit_duration_hours": 2.5,
            "fit_depth_ppt": 3.2,
            "maximum_multiple_event_statistic": 18.0,
        }
        dvt_row = {
            "sector": 54,
            "tic_id": 390988385,
            "dvt_period_days": 2.91501,
            "dvt_epoch_btjd": 3000.5001,
            "dvt_duration_hours": 2.501,
            "dvt_depth_ppt": 3.201,
            "dvt_maximum_multiple_event_statistic": 18.001,
        }
        self.assertTrue(dvt_matches_xml(xml_row, dvt_row))
        dvt_row["dvt_period_days"] = 2.92
        self.assertFalse(dvt_matches_xml(xml_row, dvt_row))

    def test_dr122_target_row_decodes_actual_contributing_sectors(self) -> None:
        vector = ["0"] * 73
        vector[54 - 14] = "1"
        vector[81 - 14] = "1"
        text = (
            "# minSector: 14 - first sector\n"
            "# maxSector: 86 - last sector\n"
            f"0000000390988385 {''.join(vector)} 1 1\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "target_info.txt"
            path.write_text(text, encoding="utf-8")
            row = parse_dr122_target_info(path, 390988385)

        self.assertEqual(row["target_info_contributing_sectors"], "54;81")
        self.assertTrue(row["target_info_made_tce"])
        self.assertTrue(row["target_info_completed_dv"])

    def test_exomast_parsers_keep_header_and_table_checks_separate(self) -> None:
        info_payload = {
            "DV Data Header": {
                "TICID": 390988385,
                "SECTORS": "s0014-s0086",
                "TPERIOD": 2.915,
                "TEPOCH": 2770.2,
                "TDEPTH": 3200.0,
                "TDUR": 2.7,
                "TSNR": 34.0,
                "MAXMES": 27.0,
                "NTRANS": 17,
                "TIMEDEL": 10.0 / 1440.0,
                "TSTART": 2769.8,
                "TSTOP": 3533.2,
            }
        }
        table_payload = {
            "data": [
                {
                    "TIME": 2769.8,
                    "LC_DETREND": 0.001,
                    "MODEL_INIT": -0.0032,
                    "SECTORS": "s0014-s0086",
                },
                {
                    "TIME": 3533.2,
                    "LC_DETREND": None,
                    "MODEL_INIT": 0.0,
                    "SECTORS": "s0014-s0086",
                },
            ]
        }
        with tempfile.TemporaryDirectory() as directory:
            info_path = Path(directory) / "info.json"
            table_path = Path(directory) / "table.json"
            info_path.write_text(json.dumps(info_payload), encoding="utf-8")
            table_path.write_text(json.dumps(table_payload), encoding="utf-8")
            info = parse_exomast_info(info_path)
            table = parse_exomast_table(table_path)

        self.assertEqual(info["exomast_scope"], "s0014-s0086")
        self.assertAlmostEqual(info["exomast_cadence_minutes"], 10.0)
        self.assertEqual(table["exomast_table_rows"], 2)
        self.assertEqual(table["exomast_table_finite_detrended_rows"], 1)
        self.assertAlmostEqual(table["exomast_table_model_depth_ppt"], 3.2)

    def test_combined_dvt_header_has_sector_vector_without_sector_keyword(
        self,
    ) -> None:
        vector = ["0"] * 82
        vector[54] = "1"
        vector[81] = "1"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "combined_dvt.fits"
            primary = fits.PrimaryHDU()
            primary.header["TICID"] = 390988385
            primary.header["DATA_REL"] = 122
            primary.header["SECTORS"] = "".join(vector)
            column = fits.Column(
                name="LC_DETREND",
                format="E",
                array=np.asarray([0.0, np.nan], dtype=np.float32),
            )
            tce = fits.BinTableHDU.from_columns([column], name="TCE_1")
            tce.header["TIMEDEL"] = 10.0 / 1440.0
            tce.header["TPERIOD"] = 2.915
            tce.header["TEPOCH"] = 2770.2
            tce.header["TDUR"] = 2.7
            tce.header["TDEPTH"] = 3200.0
            tce.header["MAXMES"] = 27.0
            fits.HDUList([primary, tce]).writeto(path)
            row = parse_dvt_header(path)

        self.assertIsNone(row["sector"])
        self.assertEqual(decode_spoc_sector_vector("".join(vector)), (54, 81))
        self.assertEqual(row["dvt_sectors_used"], "54;81")
        self.assertEqual(row["dvt_finite_detrended_rows"], 1)
        self.assertAlmostEqual(row["dvt_cadence_minutes"], 10.0)


if __name__ == "__main__":
    unittest.main()
