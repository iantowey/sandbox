package com.lucidlogic.horsetracker.presenter;

import com.lucidlogic.horsetracker.R;
import com.lucidlogic.horsetracker.model.Race;
import com.lucidlogic.horsetracker.model.Racecard;
import com.lucidlogic.horsetracker.model.dto.RaceDTO;
import com.lucidlogic.horsetracker.presenter.impl.RacePresenterImpl;
import com.lucidlogic.horsetracker.presenter.impl.RacecardPresenterImpl;
import com.lucidlogic.horsetracker.service.ParseService;
import com.lucidlogic.horsetracker.service.ParseTest;
import com.lucidlogic.horsetracker.utils.BeanTransformers;
import com.lucidlogic.horsetracker.view.RaceView;
import com.lucidlogic.horsetracker.view.RacecardView;
import com.parse.ParseException;
import com.parse.ParseQuery;

import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.mockito.Mockito;
import org.powermock.api.mockito.PowerMockito;
import org.powermock.core.classloader.annotations.PrepareForTest;
import org.powermock.modules.junit4.PowerMockRunner;
import org.robolectric.Robolectric;

import rx.Scheduler;
import rx.android.plugins.RxAndroidSchedulersHook;
import rx.android.schedulers.AndroidSchedulers;
import rx.parse.ParseObservable;
import rx.schedulers.Schedulers;
import timber.log.Timber;

import static org.mockito.Matchers.any;
import static org.powermock.api.mockito.PowerMockito.mockStatic;

/**
 * Created by ian on 01/08/16.
 */
//@PrepareForTest({AndroidSchedulers.class})
public class RacePresenterTest extends ParseTest {

    private RacePresenter racePresenter;

//    @Before
//    public void before(){
//
//        PowerMockito.stub(PowerMockito.method(AndroidSchedulers.class, "mainThread")).toReturn(Schedulers.immediate());
//        PowerMockito.stub(PowerMockito.method(Schedulers.class, "io")).toReturn(Schedulers.immediate());
//
//        racePresenter = new RacePresenterImpl();
//        RaceView raceView = Mockito.mock(RaceView.class);
//        Mockito.doAnswer(answer ->{
//            Timber.i("View invoked");
//            return null;})
//                .when(raceView)
//                .updateView(any(Race.class));
//        racePresenter.setView(raceView);
//    }

    @Test
    public void updateViewTest() throws ParseException {
        racePresenter.updateView("2wpnuNTeYC");

        RaceDTO  raceDTO = ParseObservable.find(
                ParseQuery.getQuery(RaceDTO.class)
                        .whereEqualTo("objectId", "2wpnuNTeYC")
                        .include("runners")
                        .include("runners.jockey")
                        .include("runners.horse")
                        .include("runners.trainer")
        ).toBlocking().first();

        raceDTO.pin();


    }
}
