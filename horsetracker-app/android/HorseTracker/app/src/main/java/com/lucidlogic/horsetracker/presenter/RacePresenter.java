package com.lucidlogic.horsetracker.presenter;

import com.lucidlogic.horsetracker.view.RaceView;

/**
 * Created by ian on 10/07/16.
 */
public interface RacePresenter {
    void setView(RaceView raceView);
    void updateView(String raceObjectId);
}
