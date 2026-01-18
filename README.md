# Public Night Helper

A Python command-line tool that generates astronomical data for all Saturdays in a given year. This tool is designed to help plan public astronomy observation nights by providing sunset times, moon rise/set times, moon illumination percentages, and nearby US holidays.

## Features

- Generates astronomical data for all Saturdays in a specified year
- Calculates sunset times
- Calculates moonrise and moonset times
- Computes moon illumination percentage
- Identifies US holidays within the Thursday-Monday window around each Saturday
- Outputs data to CSV format for easy analysis

## Location

The tool is configured for coordinates: **38°53'24.8"N, 96°00'08.7"W** using the **America/Chicago** timezone.

## Requirements

- Python 3.12
- holidays
- astropy

## Installation

1. Clone this repository
2. Install dependencies:

```bash
uv install
```

## Usage

```bash
uv run python main.py <year> [-o OUTPUT_FILE]
```

### Arguments

- `year`: The year to generate data for (must be between 1900 and 2100)
- `-o, --output`: Optional output CSV filename (default: `saturdays_YYYY.csv`)

### Examples

Generate data for 2024:
```bash
python main.py 2024
```

Generate data for 2025 with custom output file:
```bash
python main.py 2025 -o astronomy_2025.csv
```

## Output Format

The tool generates a CSV file with the following columns:

- `date`: Saturday date in YYYY-MM-DD format
- `sun_set_time`: Sunset time in HH:MM format (local time)
- `moon_rise_time`: Moonrise time in HH:MM format or "N/A"
- `moon_set_time`: Moonset time in HH:MM format or "N/A"
- `moon_illumination_percent`: Moon illumination as percentage (0.0-100.0)
- `holiday`: US holidays occurring Friday through Monday around the Saturday (semicolon-separated if multiple)

## How It Works

The tool uses the [Astropy](https://www.astropy.org/) library to perform accurate astronomical calculations:

1. Identifies all Saturdays in the specified year
2. For each Saturday:
   - Calculates sunset time by detecting when the sun crosses the horizon
   - Calculates moonrise and moonset times
   - Computes moon illumination based on the sun-moon elongation angle
   - Checks for US holidays within a 5-day window (Thursday-Monday)
3. Exports all data to a CSV file

## Notes

- Times are calculated for the America/Chicago timezone
- Moon illumination is calculated at sunset time (or noon if sunset is unavailable)
- The horizon crossing detection uses a 15-minute step with refinement to the nearest minute
- Holiday categories include PUBLIC, UNOFFICIAL, and GOVERNMENT holidays
