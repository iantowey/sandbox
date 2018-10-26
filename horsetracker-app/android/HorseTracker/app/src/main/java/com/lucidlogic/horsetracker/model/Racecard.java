package com.lucidlogic.horsetracker.model;

import android.databinding.BaseObservable;

import com.google.repacked.kotlin.collections.EmptyList;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;

/**
 * Created by ian on 06/07/16.
 */
public class Racecard extends BaseObservable {
    private String id;
    private List<Meeting> meetings = new ArrayList<>();
    private Date date;

    public void setId(String id) {
        this.id = id;
    }

    public List<Meeting> getMeetings() {
        return meetings;
    }

    public void setDate(Date date) {
        this.date = date;
    }

    public String getId() {
        return id;
    }

    public void setMeetings(List<Meeting> meetings) {
        this.meetings = meetings;
    }

    public Date getDate() {
        return date;
    }
}
