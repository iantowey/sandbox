package com.lucidlogic.horsetracker.model.dto;

import com.lucidlogic.horsetracker.config.Constants;
import com.parse.ParseClassName;
import com.parse.ParseObject;

import java.util.List;

/**
 * Created by ian on 06/07/16.
 */
@ParseClassName(Constants.MEETING)
public class MeetingDTO extends ParseObject{

    public String getTrack(){
        return getString(Constants.TRACK);
    }
    public String getTrackSurface(){
        return getString(Constants.TRACK_SURFACE);
    }
    public String getTrackGoing(){
        return getString(Constants.TRACK_GOING);
    }
    public List<RaceDTO> getRaces(){
        return getList(Constants.RACES);
    }

    public void getTrack(String track){
        put(Constants.TRACK, track);
    }
    public void getTrackSurface(String trackSurface){
        put(Constants.TRACK_SURFACE, trackSurface);
    }
    public void getTrackGoing(String trackGoing){
        put(Constants.TRACK_GOING, trackGoing);
    }
    public void getRaces(List<RaceDTO> races){
        put(Constants.RACES, races);
    }

}
