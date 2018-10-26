package com.lucidlogic.horsetracker.model;

import android.databinding.BaseObservable;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by ian on 06/07/16.
 */
public class Meeting extends BaseObservable {
    private String going;
    private String id;
    private String surface;
    private String track;
    private List<Race> races = new ArrayList<>();

    public void setGoing(String going) {
        this.going = going;
    }

    public void setId(String id) {
        this.id = id;
    }

    public void setSurface(String surface) {
        this.surface = surface;
    }

    public void setTrack(String track) {
        this.track = track;
    }


    public List<Race> getRaces() {
        return races;
    }

    public void setRaces(List<Race> races) {
        this.races = races;
    }

    public String getTrack() {
        return track;
    }

    public String getSurface() {
        return surface;
    }

    public String getId() {
        return id;
    }

    public String getGoing() {
        return going;
    }
}
