package com.lucidlogic.horsetracker.presenter.impl;

import com.lucidlogic.horsetracker.presenter.RacePresenter;
import com.lucidlogic.horsetracker.presenter.RacecardPresenter;
import com.lucidlogic.horsetracker.service.ParseService;
import com.lucidlogic.horsetracker.utils.BeanTransformers;
import com.lucidlogic.horsetracker.view.RaceView;
import com.lucidlogic.horsetracker.view.RacecardView;
import com.parse.ParseException;

import java.util.Date;

import rx.android.schedulers.AndroidSchedulers;
import rx.schedulers.Schedulers;
import timber.log.Timber;

/**
 * Created by ian on 10/07/16.
 */
public class RacePresenterImpl implements RacePresenter {

    RaceView raceView;

    @Override
    public void setView(RaceView raceView) {
        this.raceView = raceView;
    }

    @Override
    public void updateView(String raceObjectId) {
        ParseService.getRace(raceObjectId)
                .observeOn(AndroidSchedulers.mainThread())
                .subscribeOn(Schedulers.io())
                .map(raceDTO -> {
//                    raceDTO.pinInBackground();
                    return BeanTransformers.raceFromRaceDTO(raceDTO, true);
                })
                .subscribe(race -> {
                            Timber.i("Item retrieved");
                            raceView.updateView(race);
                        },
                        e -> {
                            Timber.e(e,e.getMessage());
                            e.printStackTrace();
                        },
                        () -> {
                            Timber.i("Item retrival complete");
                        });

    }
}
