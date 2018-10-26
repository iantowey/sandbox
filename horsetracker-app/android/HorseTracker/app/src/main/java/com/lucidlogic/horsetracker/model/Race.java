package com.lucidlogic.horsetracker.model;

import android.databinding.BaseObservable;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by ian on 06/07/16.
 */
public class Race extends BaseObservable {

    private String id;
    private String abandoned;
    private String name;
    private String time;
    private List<Runner> runners = new ArrayList<>();

    public void setAbandoned(String abandoned) {
        this.abandoned = abandoned;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setTime(String time) {
        this.time = time;
    }

    public String getAbandoned() {
        return abandoned;
    }

    public String getName() {
        return name;
    }

    public String getTime() {
        return time;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public List<Runner> getRunners() {
        return runners;
    }

    public void setRunners(List<Runner> runners) {
        this.runners = runners;
    }
}
