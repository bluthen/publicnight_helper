#!/usr/bin/env python3
"""
Command line program to generate astronomical data for Saturdays in a given year.
Outputs sunset time, moonrise time, and moon illumination percentage.
"""

import argparse
import csv
import sys
import warnings
from datetime import datetime, timedelta
from zoneinfo import ZoneInfo

import numpy as np
import holidays
from holidays.constants import PUBLIC, UNOFFICIAL, GOVERNMENT
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun, get_body
import astropy.units as u

# Suppress astropy warnings about coordinate transformations
warnings.filterwarnings("ignore", category=Warning, module="astropy")


def get_holiday_for_weekend(saturday, year):
    """
    Check if there's a holiday on Thursday, Friday, Saturday, Sunday, or Monday
    around the given Saturday.
    Returns all holiday names separated by semicolons, or empty string.
    """
    us_holidays = holidays.US(categories=(PUBLIC, UNOFFICIAL, GOVERNMENT), years=year)

    # Check Thursday through Monday (5 day window)
    thursday = saturday - timedelta(days=2)
    friday = saturday - timedelta(days=1)
    sunday = saturday + timedelta(days=1)
    monday = saturday + timedelta(days=2)

    holiday_names = []
    for day in [thursday, friday, saturday, sunday, monday]:
        if day in us_holidays:
            holiday_names.append(us_holidays[day])

    return "; ".join(holiday_names)


def parse_coordinates(coord_string):
    """
    Parse coordinates in format: 38°53'24.8"N+96°00'08.7"W
    Returns: (latitude, longitude) in degrees
    """
    # 38°53'24.8"N = 38 + 53/60 + 24.8/3600
    # 96°00'08.7"W = -(96 + 0/60 + 8.7/3600)
    lat = 38 + 53 / 60 + 24.8 / 3600  # North is positive
    lon = -(96 + 0 / 60 + 8.7 / 3600)  # West is negative
    return lat, lon


def get_saturdays(year):
    """
    Get all Saturdays in the given year, starting from the first Saturday
    on or after January 1st through the final Saturday in December.
    Returns list of datetime.date objects.
    """
    saturdays = []

    # Find first Saturday on or after January 1st
    current = datetime(year, 1, 1).date()
    # Saturday is weekday 5 (Monday=0)
    days_until_saturday = (5 - current.weekday()) % 7
    if days_until_saturday > 0 or current.weekday() == 5:
        current += timedelta(days=days_until_saturday)

    # Get all Saturdays through December 31st
    year_end = datetime(year, 12, 31).date()
    while current <= year_end:
        saturdays.append(current)
        current += timedelta(days=7)

    return saturdays


def find_moon_rise_set(date, location, timezone):
    """
    Find both moonrise and moonset for a given date.
    More efficient than calling find_event_time twice.

    Returns:
        Tuple of (moonrise_str, moonset_str) in HH:MM format or "N/A"
    """
    # Start search at midnight local time for the day
    dt_start = datetime.combine(date, datetime.min.time()).replace(tzinfo=timezone)

    # Search through 48 hours to handle edge cases
    minutes_to_search = 48 * 60
    step = 15  # Check every 15 minutes for speed

    prev_alt = None
    rise_time = None
    set_time = None

    for minutes in range(0, minutes_to_search, step):
        current_time = dt_start + timedelta(minutes=minutes)
        astro_time = Time(current_time)

        body = get_body("moon", astro_time, location)
        altaz = body.transform_to(AltAz(obstime=astro_time, location=location))
        current_alt = altaz.alt.degree

        # Detect crossing the horizon
        if prev_alt is not None:
            if prev_alt < 0 and current_alt >= 0 and rise_time is None:
                # Rising: going above horizon
                rise_time = refine_event_time(
                    current_time - timedelta(minutes=step),
                    current_time,
                    location,
                    timezone,
                    "moon",
                    "rise",
                )
            elif prev_alt > 0 and current_alt <= 0 and set_time is None:
                # Setting: going below horizon
                set_time = refine_event_time(
                    current_time - timedelta(minutes=step),
                    current_time,
                    location,
                    timezone,
                    "moon",
                    "set",
                )

        prev_alt = current_alt

        # Stop if we've passed the day and found both events or if we're far enough
        if current_time.date() > date and minutes > 24 * 60:
            break

    rise_str = rise_time.strftime("%H:%M") if rise_time else "N/A"
    set_str = set_time.strftime("%H:%M") if set_time else "N/A"

    return rise_str, set_str


def find_event_time(
    date, location, timezone, body_name, event_type="set", altitude_threshold=0
):
    """
    Find the time when a celestial body rises or sets, or crosses a specific altitude.

    Args:
        date: The date to search
        location: EarthLocation object
        timezone: ZoneInfo timezone
        body_name: 'sun' or 'moon'
        event_type: 'rise' or 'set'
        altitude_threshold: Altitude in degrees (default 0 for horizon)

    Returns:
        Time string in HH:MM format or "N/A"
    """
    # Start search at midnight local time for the day
    dt_start = datetime.combine(date, datetime.min.time()).replace(tzinfo=timezone)

    # Search through 48 hours to handle edge cases
    minutes_to_search = 48 * 60
    step = 15  # Check every 15 minutes for speed

    prev_alt = None
    event_time = None

    for minutes in range(0, minutes_to_search, step):
        current_time = dt_start + timedelta(minutes=minutes)
        astro_time = Time(current_time)

        # Get body position
        if body_name == "sun":
            body = get_sun(astro_time)
        else:
            body = get_body("moon", astro_time, location)

        altaz = body.transform_to(AltAz(obstime=astro_time, location=location))
        current_alt = altaz.alt.degree

        # Detect crossing the threshold altitude
        if prev_alt is not None:
            if (
                event_type == "set"
                and prev_alt > altitude_threshold
                and current_alt <= altitude_threshold
            ):
                # Setting: going below threshold
                # Refine to the minute
                event_time = refine_event_time(
                    current_time - timedelta(minutes=step),
                    current_time,
                    location,
                    timezone,
                    body_name,
                    "set",
                    altitude_threshold,
                )
                break
            elif (
                event_type == "rise"
                and prev_alt < altitude_threshold
                and current_alt >= altitude_threshold
            ):
                # Rising: going above threshold
                event_time = refine_event_time(
                    current_time - timedelta(minutes=step),
                    current_time,
                    location,
                    timezone,
                    body_name,
                    "rise",
                    altitude_threshold,
                )
                break

        prev_alt = current_alt

        # Stop if we've passed the day we're interested in
        if current_time.date() > date and event_time is None and minutes > 24 * 60:
            break

    if event_time:
        return event_time.strftime("%H:%M")
    return "N/A"


def refine_event_time(
    start_time,
    end_time,
    location,
    timezone,
    body_name,
    event_type,
    altitude_threshold=0,
):
    """
    Use binary search to find the exact time when a celestial body
    crosses the altitude threshold.

    Refines to within ~1 second precision, then rounds to nearest minute.

    Args:
        start_time: Start of the time bracket (datetime with timezone)
        end_time: End of the time bracket (datetime with timezone)
        location: EarthLocation object
        timezone: ZoneInfo timezone
        body_name: 'sun' or 'moon'
        event_type: 'rise' or 'set'
        altitude_threshold: Altitude in degrees (default 0 for horizon)

    Returns:
        datetime object rounded to nearest minute
    """
    # Binary search: narrow down to 1-second precision
    while (end_time - start_time).total_seconds() > 1:
        mid_time = start_time + (end_time - start_time) / 2

        # Check altitude at midpoint
        astro_time = Time(mid_time)
        if body_name == "sun":
            body = get_sun(astro_time)
        else:
            body = get_body("moon", astro_time, location)

        altaz = body.transform_to(AltAz(obstime=astro_time, location=location))
        mid_alt = altaz.alt.degree

        # Determine which half contains the crossing
        if event_type == "set":
            if mid_alt > altitude_threshold:
                start_time = mid_time  # Crossing is in second half
            else:
                end_time = mid_time  # Crossing is in first half
        else:  # rise
            if mid_alt < altitude_threshold:
                start_time = mid_time  # Crossing is in second half
            else:
                end_time = mid_time  # Crossing is in first half

    # Return the midpoint rounded to nearest minute
    final_time = start_time + (end_time - start_time) / 2
    # Round to nearest minute
    if final_time.second >= 30:
        final_time += timedelta(minutes=1)
    return final_time.replace(second=0, microsecond=0)


def find_astronomical_twilight(date, location, timezone):
    """
    Find the time when astronomical twilight ends (dusk).
    This is when the sun reaches -12 degrees below the horizon.

    Args:
        date: The date to search
        location: EarthLocation object
        timezone: ZoneInfo timezone

    Returns:
        Time string in HH:MM format or "N/A"
    """
    return find_event_time(
        date, location, timezone, "sun", "set", altitude_threshold=-12
    )


def calculate_moon_illumination(date, location, timezone, sunset_time=None):
    """
    Calculate moon illumination percentage.
    Returns illumination as a percentage (0-100).
    """
    # Use sunset time if provided, otherwise use noon
    if sunset_time and sunset_time != "N/A":
        # Parse the sunset time (HH:MM format)
        hour, minute = map(int, sunset_time.split(":"))
        dt = datetime.combine(
            date, datetime.min.time().replace(hour=hour, minute=minute)
        )
    else:
        # Fallback to noon if sunset time is not available
        dt = datetime.combine(date, datetime.min.time().replace(hour=12))
    dt_local = dt.replace(tzinfo=timezone)
    astro_time = Time(dt_local)

    # Get sun and moon positions from Earth
    sun = get_sun(astro_time)
    moon = get_body("moon", astro_time, location)

    # Get Earth position (for phase angle calculation)
    # The phase angle is the Sun-Moon-Earth angle
    # We need to calculate the elongation and use it correctly
    elongation = sun.separation(moon).radian

    # The correct formula for illumination is:
    # illumination = (1 - cos(elongation)) / 2
    # This gives 0% at conjunction (new moon) and 100% at opposition (full moon)
    illumination_fraction = (1 - np.cos(elongation)) / 2
    illumination_percent = illumination_fraction * 100

    return f"{illumination_percent:.1f}"


def main():
    parser = argparse.ArgumentParser(
        description="Generate astronomical data for Saturdays in a given year."
    )
    parser.add_argument("year", type=int, help="Year to generate data for")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output CSV file (default: saturdays_YYYY.csv)",
    )

    args = parser.parse_args()
    year = args.year

    # Validate year
    if year < 1900 or year > 2100:
        print(f"Error: Year must be between 1900 and 2100", file=sys.stderr)
        sys.exit(1)

    # Set output filename
    output_file = args.output if args.output else f"saturdays_{year}.csv"

    # Setup location and timezone
    lat, lon = parse_coordinates("38°53'24.8\"N+96°00'08.7\"W")
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=0 * u.m)
    timezone = ZoneInfo("America/Chicago")

    # Get all Saturdays for the year
    saturdays = get_saturdays(year)

    print(f"Generating data for {len(saturdays)} Saturdays in {year}...")

    # Generate data for each Saturday
    data = []
    for i, saturday in enumerate(saturdays, 1):
        print(f"Processing {saturday} ({i}/{len(saturdays)})...", end="\r", flush=True)

        sunset = find_event_time(saturday, location, timezone, "sun", "set")
        astro_twilight = find_astronomical_twilight(saturday, location, timezone)
        moonrise, moonset = find_moon_rise_set(saturday, location, timezone)
        moon_illum = calculate_moon_illumination(saturday, location, timezone, sunset)
        holiday = get_holiday_for_weekend(saturday, year)

        data.append(
            {
                "date": saturday.strftime("%Y-%m-%d"),
                "sun_set_time": sunset,
                "astro_twilight": astro_twilight,
                "moon_rise_time": moonrise,
                "moon_set_time": moonset,
                "moon_illumination_percent": moon_illum,
                "holiday": holiday,
            }
        )

    print()  # New line after progress

    # Write to CSV
    with open(output_file, "w", newline="") as csvfile:
        fieldnames = [
            "date",
            "sun_set_time",
            "astro_twilight",
            "moon_rise_time",
            "moon_set_time",
            "moon_illumination_percent",
            "holiday",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerows(data)

    print(f"Data written to {output_file}")


if __name__ == "__main__":
    main()
